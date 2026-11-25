"""InSight SEIS waveform preprocessing and guarded OpenSees excitation utilities.

Units: seconds, meters, kilograms, newtons; acceleration m/s^2.
A surface-observed SEIS record is NOT a verified outcrop/base input motion.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np


INSIGHT_ACC_PATH = Path(__file__).resolve().parents[2] / "structure" / "MARS_DATA" / "S1222a_ACC.mseed"


def inspect_mseed_headers(paths):
    """Read MiniSEED v2 fixed headers without decoding samples or requiring ObsPy."""
    import struct
    from datetime import datetime, timedelta, timezone
    results = []
    for filename in paths:
        path = Path(filename)
        with path.open('rb') as f:
            b = f.read(64)
        if len(b) < 56 or b[6:7] not in (b'D', b'R', b'Q', b'M'):
            raise ValueError(f'MiniSEED v2 header is invalid: {path}')
        year, day, hour, minute, second, unused, tenth_ms = struct.unpack('>HHBBBBH', b[20:30])
        sample_count, rate_factor, rate_multiplier = struct.unpack('>Hhh', b[30:36])
        if rate_factor == 0 or rate_multiplier == 0:
            sample_rate = None
        else:
            f = rate_factor if rate_factor > 0 else 1 / abs(rate_factor)
            m = rate_multiplier if rate_multiplier > 0 else 1 / abs(rate_multiplier)
            sample_rate = f * m
        dt = datetime(year,1,1,tzinfo=timezone.utc) + timedelta(days=day-1,hours=hour,minutes=minute,seconds=second,microseconds=tenth_ms*100)
        results.append({'file':str(path),'bytes':path.stat().st_size,
                        'network':b[18:20].decode('ascii').strip(),
                        'station':b[8:13].decode('ascii').strip(),
                        'location':b[13:15].decode('ascii').strip(),
                        'channel':b[15:18].decode('ascii').strip(),
                        'first_record_start_utc':dt.isoformat(),
                        'first_record_samples':sample_count,
                        'sample_rate_hz':sample_rate,
                        'note':'Header-only inspection; data samples and event alignment not yet verified.'})
    channels = {item['channel'] for item in results}
    if channels != {'BHU','BHV','BHW'}:
        raise ValueError(f'Expected BHU/BHV/BHW, got {sorted(channels)}')
    return results


def preprocess_mseed(mseed_paths, stationxml, output_csv, freqmin=0.05, freqmax=2.0, target_rate=10.0):
    """Correct instrument response to acceleration, rotate sensor U/V/W -> Z/N/E."""
    try:
        from obspy import read, read_inventory
        from obspy.signal.rotate import rotate2zne
    except ImportError as exc:
        raise RuntimeError('ObsPy 필요: pip install obspy') from exc
    if len(mseed_paths) != 3:
        raise ValueError('같은 이벤트, 동일 대역의 U/V/W 세 파일이 필요합니다.')
    inventory = read_inventory(str(stationxml))
    stream = None
    for filename in mseed_paths:
        s = read(str(filename))
        stream = s if stream is None else stream + s
    stream.merge(method=1, fill_value=None)
    if len(stream) != 3:
        raise ValueError(f'병합된 채널 수가 3이 아닙니다: {len(stream)}. 누락/중복을 확인하세요.')
    ids = [tr.id for tr in stream]
    suffixes = [tr.stats.channel[-1].upper() for tr in stream]
    if set(suffixes) != {'U','V','W'}:
        raise ValueError(f'U/V/W 세 성분이 필요합니다: {ids}')
    starts = [tr.stats.starttime for tr in stream]
    ends = [tr.stats.endtime for tr in stream]
    start, end = max(starts), min(ends)
    if end <= start:
        raise ValueError('세 성분의 공통 관측 구간이 없습니다.')
    if not (0 < freqmin < freqmax < target_rate / 2):
        raise ValueError('0 < freqmin < freqmax < Nyquist 조건이 필요합니다.')
    stream.trim(starttime=start, endtime=end)
    azdip = {}
    for tr in stream:
        selection = inventory.select(network=tr.stats.network, station=tr.stats.station,
                                     location=tr.stats.location, channel=tr.stats.channel, time=start)
        try:
            channel = selection[0][0][0]
            azdip[tr.id] = (float(channel.azimuth), float(channel.dip))
        except (IndexError, AttributeError, TypeError) as exc:
            raise RuntimeError(f'채널 방향 메타데이터 부족: {tr.id}') from exc
        tr.detrend('linear')
        tr.taper(max_percentage=0.05)
        tr.remove_response(inventory=inventory, output='ACC',
                           pre_filt=(freqmin*0.5, freqmin, freqmax, min(freqmax*1.5,target_rate/2*0.95)),
                           water_level=60, zero_mean=False, taper=False)
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
        tr.resample(target_rate, no_filter=False)
    stream.trim(starttime=max(tr.stats.starttime for tr in stream),
                endtime=min(tr.stats.endtime for tr in stream))
    count = min(len(tr.data) for tr in stream)
    if count < 3:
        raise RuntimeError('유효한 공통 샘플 부족')
    comps = []
    for tr in stream:
        if np.ma.isMaskedArray(tr.data) or not np.all(np.isfinite(tr.data[:count])):
            raise ValueError(f'결측/비정상 샘플: {tr.id}')
        comps.append(tr.data[:count].astype(float))
    z, n, e = rotate2zne(comps[0], *azdip[ids[0]], comps[1], *azdip[ids[1]], comps[2], *azdip[ids[2]])
    time = np.arange(count, dtype=float) / float(target_rate)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_csv, np.column_stack((time,e,n,z)), delimiter=',',
               header='time_s,east_m_s2,north_m_s2,up_m_s2', comments='')
    metadata = {'data_type':'surface_observed_acceleration_NOT_base_input',
                'units':'m/s^2', 'coordinate_system':'E,N,Up', 'dt_s':1/target_rate,
                'samples':count, 'start_utc':str(start), 'stationxml':str(stationxml),
                'channels':ids,'frequency_band_hz':[freqmin,freqmax],
                'caution':'Surface SEIS data must not be treated as verified bedrock input motion.'}
    output_csv.with_suffix('.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False),encoding='utf-8')
    return output_csv



def preprocess_acc_mseed(filename, output_csv, freqmin=0.05, freqmax=2.0, target_rate=10.0):
    """Read the supplied response-corrected S1222a ACC (m/s²), already Z/N/E.

    MiniSEED does not encode physical units: this entry point is only for the
    supplied calibrated ACC product, never raw counts. Do not remove response
    again or rotate these geographic channels with the old UVW StationXML.
    """
    from obspy import read
    stream = read(str(filename))
    stream.merge(method=1, fill_value=None)
    if len(stream) != 3 or {tr.stats.channel for tr in stream} != {'BHE', 'BHN', 'BHZ'}:
        raise ValueError('ACC 파일에는 BHE/BHN/BHZ 세 성분이 필요합니다.')
    if len({(tr.stats.network, tr.stats.station, tr.stats.location) for tr in stream}) != 1:
        raise ValueError('ACC 성분의 관측소가 다릅니다.')
    ref = stream[0].stats.copy()
    if not (0 < freqmin < freqmax < min(target_rate, ref.sampling_rate) / 2):
        raise ValueError('주파수 범위는 원본 및 출력 Nyquist보다 작아야 합니다.')
    for tr in stream:
        if (tr.stats.starttime != ref.starttime or tr.stats.npts != ref.npts
                or tr.stats.sampling_rate != ref.sampling_rate):
            raise ValueError('ACC 성분의 시간축이 일치하지 않습니다.')
        if np.ma.isMaskedArray(tr.data) or not np.all(np.isfinite(tr.data)):
            raise ValueError(f'결측/비정상 ACC 샘플: {tr.id}')
        tr.detrend('linear')
        tr.taper(max_percentage=0.05)
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
        if tr.stats.sampling_rate != target_rate:
            tr.resample(target_rate, no_filter=False)
    xyz = np.column_stack([stream.select(channel=c)[0].data for c in ('BHE', 'BHN', 'BHZ')])
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_csv, np.column_stack((np.arange(len(xyz))/target_rate, xyz)),
               delimiter=',', header='time_s,east_m_s2,north_m_s2,up_m_s2', comments='')
    metadata = dict(data_type='surface_observed_acceleration_NOT_base_input',
                    units='m/s^2', coordinate_system='E,N,Up', dt_s=1/target_rate,
                    samples=len(xyz), start_utc=str(ref.starttime),
                    source_file=str(Path(filename).resolve()), event='S1222a',
                    channels=[stream.select(channel=c)[0].id for c in ('BHE','BHN','BHZ')],
                    frequency_band_hz=[freqmin, freqmax],
                    response_removed_in_this_step=False,
                    units_basis='Supplied S1222a_ACC product assumed calibrated in m/s^2; MiniSEED has no unit metadata',
                    caution='Surface SEIS data must not be treated as verified bedrock input motion.')
    output_csv.with_suffix('.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding='utf-8')
    return output_csv


def load_acceleration_csv(filename):
    path=Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f'가속도 CSV 없음: {path}')
    data=np.genfromtxt(path,delimiter=',',names=True,dtype=float,encoding='utf-8')
    if data.ndim == 0:
        raise ValueError('샘플 부족')
    expected=('time_s','east_m_s2','north_m_s2','up_m_s2')
    if tuple(data.dtype.names or ()) != expected:
        raise ValueError(f'CSV 컬럼은 {expected} 여야 합니다.')
    values=np.column_stack([data[k] for k in expected])
    if len(values)<3 or not np.all(np.isfinite(values)):
        raise ValueError('가속도 데이터에 결측값이 있거나 샘플이 부족합니다.')
    intervals=np.diff(values[:,0]); dt=float(np.median(intervals))
    if dt<=0 or not np.allclose(intervals,dt,rtol=1e-3,atol=1e-7):
        raise ValueError('시간은 일정한 간격으로 증가해야 합니다.')
    if abs(values[0,0])>dt*1e-3:
        raise ValueError('시간 열은 0초에서 시작해야 합니다.')
    return dt,values[:,1:]


def run_uniform_excitation(model_data, csv_path, *, acknowledge_surface_as_base=False,
                           scale=1.0, max_steps=None, output_csv=None,
                           animation_html=None, animation_stride=10, visual_scale=100.0):
    """EXPERIMENTAL uniform input; not a validated Mars site-response simulation."""
    if not acknowledge_surface_as_base:
        raise ValueError('SEIS 지표 관측파형을 기반암 운동으로 취급하는 가정을 명시적으로 허용해야 합니다.')
    if model_data.get('gravity_summary') is None:
        raise ValueError('동적해석 전 중력 초기화가 필요합니다.')
    if not np.isfinite(scale) or scale<=0:
        raise ValueError('scale은 양수여야 합니다.')
    if (animation_stride is not None and int(animation_stride) < 1) or (max_steps is not None and int(max_steps) < 1):
        raise ValueError('animation_stride and max_steps must be >= 1')
    from openseespy import opensees as ops
    dt, xyz=load_acceleration_csv(csv_path)
    sidecar = Path(csv_path).with_suffix('.json')
    if not sidecar.exists():
        raise ValueError('가속도 CSV의 처리 메타데이터 JSON이 없습니다. 검증되지 않은 수치 사용을 차단합니다.')
    meta = json.loads(sidecar.read_text(encoding='utf-8'))
    if meta.get('data_type') != 'surface_observed_acceleration_NOT_base_input' or meta.get('coordinate_system') != 'E,N,Up' or meta.get('units') != 'm/s^2':
        raise ValueError('CSV 출처/단위/좌표계가 기기보정된 ENU 지표 기록과 일치하지 않습니다.')
    model_data['motion_metadata'] = meta
    if model_data.get('boundary_condition_mode') != 'experimental_fixed_base_free_sides':
        raise ValueError('측면 롤러가 수평 동적 자유도를 제한합니다. 실험용 측면 자유 경계로 재생성해야 합니다.')
    # XYZ corresponds to East/North/Up ONLY when ground FEM x/y are aligned ENU.
    # Horizontal axes may differ in rotated/uncalibrated DEM. Stop by default.
    if not model_data.get('dynamic_axes_verified_enu', False):
        raise ValueError('지반 FEM의 X/Y 방향이 East/North와 정렬되었는지 검증되지 않았습니다.')
    ops.wipeAnalysis()
    for i in range(3):
        tag=5001+i
        # Use in-memory values; OpenSees Path accepts list via -values.
        ops.timeSeries('Path',tag,'-dt',dt,'-values',*map(float,xyz[:,i]),'-factor',float(scale))
        ops.pattern('UniformExcitation',6001+i,i+1,'-accel',tag)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    try: ops.system('UmfPack')
    except Exception: ops.system('BandGeneral')
    ops.test('NormDispIncr',1e-7,60)
    ops.algorithm('Newton')
    ops.integrator('Newmark',0.5,0.25)
    ops.analysis('Transient')
    surface_nodes=model_data['surface_nodes']
    center=np.mean([model_data['node_coords'][t][:2] for t in surface_nodes],axis=0)
    node=min(surface_nodes,key=lambda t:np.linalg.norm(np.array(model_data['node_coords'][t][:2])-center))
    max_n=len(xyz)-1 if max_steps is None else min(int(max_steps),len(xyz)-1)
    if animation_stride is None:
        animation_stride = max(1, int(np.ceil(max_n/180)))
    history=[]
    snapshots=[]
    ordered_tags=sorted(model_data['node_coords'])
    if animation_html is not None:
        snapshots.append({'time': float(ops.getTime()),
                          'disp': np.asarray([ops.nodeDisp(t) for t in ordered_tags], dtype=np.float64)})
    dynamic_started = time.perf_counter()
    for step in range(max_n):
        ok=ops.analyze(1,dt)
        if ok!=0:
            raise RuntimeError(f'동적 수렴 실패: step {step+1}/{max_n}; 결과를 성공으로 사용하면 안 됩니다.')
        history.append([ops.getTime(),*ops.nodeDisp(node),*ops.nodeAccel(node)])
        if (step+1) % max(1, max_n//10) == 0 or step == max_n-1:
            elapsed = time.perf_counter() - dynamic_started
            remaining = elapsed * (max_n-step-1)/(step+1)
            print(f"동적 해석 {step+1}/{max_n} · t={ops.getTime():.3f} s · 경과 {elapsed:.1f}s · 예상 잔여 {remaining:.1f}s", flush=True)
        if animation_html is not None and ((step+1) % int(animation_stride) == 0 or step == max_n-1):
            snapshots.append({'time': float(ops.getTime()),
                              'disp': np.asarray([ops.nodeDisp(t) for t in ordered_tags], dtype=np.float64)})
    if animation_html is not None:
        if int(animation_stride) < 1:
            raise ValueError('animation_stride must be >= 1')
        from insight_dynamic_viewer import save_dynamic_viewer
        save_dynamic_viewer(model_data,snapshots,history,node,animation_html,visual_scale)
    if output_csv:
        np.savetxt(output_csv,np.array(history),delimiter=',',comments='',
                   header='time_s,ux_m,uy_m,uz_m,ax_relative_m_s2,ay_relative_m_s2,az_relative_m_s2')
    return {'completed_steps':max_n,'dt_s':dt,'sample_node':int(node),'experimental':True}


def main():
    p=argparse.ArgumentParser(description='InSight 3-component MiniSEED -> calibrated surface acceleration CSV')
    p.add_argument('--mseed',type=Path,nargs=3,required=True,metavar='FILE')
    p.add_argument('--stationxml',type=Path,required=True)
    p.add_argument('--output',type=Path,default=Path('insight_acceleration_enu.csv'))
    p.add_argument('--freqmin',type=float,default=0.05)
    p.add_argument('--freqmax',type=float,default=2.0)
    p.add_argument('--rate',type=float,default=10.0)
    a=p.parse_args()
    print('Created:',preprocess_mseed(a.mseed,a.stationxml,a.output,a.freqmin,a.freqmax,a.rate))



def select_motion_window(csv_path, output_csv, duration=None, start=None):
    """Select a calibrated ENU interval; default centre is vector-acceleration peak.

    A 5% cosine taper at both window ends avoids an abrupt base input jump.
    Time is reset to zero and source offsets are retained in the JSON sidecar.
    """
    dt, xyz = load_acceleration_csv(csv_path)
    if duration is None:
        duration = (len(xyz)-1)*dt
        if start is None:
            start = 0.0
    if not np.isfinite(duration) or duration < 2*dt:
        raise ValueError('motion-duration must be finite and at least two time steps')
    count = min(len(xyz), int(round(duration/dt))+1)
    peak = int(np.argmax(np.linalg.norm(xyz, axis=1)))
    if start is None:
        first = max(0, min(peak-count//2, len(xyz)-count))
    else:
        if not np.isfinite(start) or start < 0 or start > (len(xyz)-3)*dt:
            raise ValueError('motion-start must be inside the record, with >=3 samples remaining')
        first = int(round(start/dt))
        count = min(count, len(xyz)-first)
    selected = xyz[first:first+count].copy()
    nedge = max(2, int(np.ceil(0.05*count)))
    ramp = 0.5*(1-np.cos(np.linspace(0, np.pi, nedge)))
    selected[:nedge] *= ramp[:,None]
    selected[-nedge:] *= ramp[::-1,None]
    meta = json.loads(Path(csv_path).with_suffix('.json').read_text(encoding='utf-8'))
    if meta.get('coordinate_system') != 'E,N,Up' or meta.get('units') != 'm/s^2':
        raise ValueError('Expected calibrated ENU acceleration metadata')
    meta.update(source_csv=str(Path(csv_path).resolve()), source_start_utc=meta.get('start_utc'),
                source_offset_s=first*dt, source_end_s=(first+count-1)*dt,
                source_peak_s=peak*dt, window_selection='peak_centered' if start is None else 'manual',
                window_taper_fraction=0.05, samples=count)
    if meta.get('start_utc'):
        from datetime import datetime, timedelta
        meta['start_utc'] = (datetime.fromisoformat(meta['start_utc'].replace('Z','+00:00'))
                             + timedelta(seconds=first*dt)).isoformat()
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_csv, np.column_stack((np.arange(count)*dt, selected)), delimiter=',',
               header='time_s,east_m_s2,north_m_s2,up_m_s2', comments='')
    output_csv.with_suffix('.json').write_text(json.dumps(meta, indent=2, ensure_ascii=False),encoding='utf-8')
    print(f'입력 구간: 원기록 {first*dt:.2f}–{(first+count-1)*dt:.2f} s '
          f'(관찰 {(count-1)*dt:.2f} s, 전체 최대 가속도 {peak*dt:.2f} s)', flush=True)
    return output_csv


if __name__ == "__main__":
    main()
