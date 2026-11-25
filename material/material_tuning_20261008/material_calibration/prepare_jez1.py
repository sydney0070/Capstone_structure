"""Published strength inputs, not stress-strain calibration or raw experimental targets."""
from pathlib import Path
from math import tan,radians
import json,sys
B=Path(__file__).resolve().parent;root=B/'mars_jez1_strength';root.mkdir(exist_ok=True)
rows=[]
for surface,peak,res,cpeak,cres in [('soil',40.6,39.7,3.1,3.2),('smooth_steel',18.9,16.7,1.5,1.5),('rough_steel',36.2,33.2,4.3,2.8)]:
 for state,angle,cohesion in [('peak',peak,cpeak),('residual',res,cres)]:
  rows.append(dict(surface=surface,state=state,angle_deg=angle,cohesion_or_adhesion_kPa=cohesion,friction_coefficient=tan(radians(angle))))
out=dict(status='LITERATURE_STRENGTH_TRANSFER_ONLY_NOT_CALIBRATED',source_doi='10.3390/app11157052',source_locations='Table 4; Figure 10 and section 4.3',license='CC-BY-4.0',material='JEZ-1 simulant tested on Earth; not in-situ Mars soil',density_tonne_m3=1.84,relative_density=.8628,normal_stresses_kPa=[25,50,75],inputs=rows,not_identified=['K','G','dilatancy','hardening','contact_Kn','contact_Kt','low-stress extrapolation'],model_caveat='DruckerPrager compression-meridian mapping is not direct-shear calibration. Existing zeroLengthContactASDimplex coupon has no adhesion, so mu alone cannot reproduce these envelopes; no targets are manufactured from the fitted envelopes.')
# 2026-10-19: distinguish soil and interface test conditions (paper section 3.2).
# Preserve the legacy field as the soil-only list; no fitted parameters are changed.
out['normal_stresses_kPa_scope']='soil-soil direct shear only; see interface_normal_stresses_kPa'
out['interface_normal_stresses_kPa']=[25,50,75,100]
out['shear_test_conditions']=dict(source='sections 2.2 and 3.2, PDF pages 4-5',shear_rate_mm_min=.5,max_horizontal_displacement_mm=10,upper_box_plan_mm=[100,100],mean_soil_sample_height_mm=31.59,mean_interface_sample_height_mm=27.58,roughness_Ra_um=dict(smooth=1.0,rough=45.1),rough_surface='2 mm epoxy/Fontainebleau sand layer; not a generic steel finish')
out['oedometer_published_summary']=dict(status='LITERATURE_SUMMARY_NOT_FITTED',source='sections 3.1 and 4.1, Figure 5, PDF pages 4 and 6',dry_sample=True,loading_effective_vertical_stresses_kPa=[12.5,25,50,100,200,400,800],unloading_endpoint_kPa=50,compression_index_Cc=.131,swelling_index_Cs=.030,caveat='Not directly K or G; curve digitization and model suitability review remain. Not a permeability or time-dependent consolidation calibration.')
(root/'published_inputs.json').write_text(json.dumps(out,indent=2)+'\n')
(root/'README.md').write_text('''# JEZ-1 지반·강판 접촉의 문헌 강도 입력

상태: **문헌 강도값 반영 단계, 독립 곡선 보정 미완료**.

[원 연구: Yin 등 (2021)](https://doi.org/10.3390/app11157052)의 Table 4와 section 4.3에서 지반의 c·φ, 강판 접촉의 부착력·마찰각을 읽어 `published_inputs.json`에 기록했다. 수치 범위·시료 밀도·응력 조건도 함께 보관했다. `prepare_jez1.py`로 재생성한다.

지반 후보는 `DruckerPrager` + `SSPbrick`이다. `materials.py`의 c·φ 변환은 삼축압축 자오선 대응이며 직접전단 시험 재현을 보증하지 않는다. K·G·팽창각·경화는 아직 정하지 않았다. peak/residual은 두 상태의 값이지 이를 연결하는 연화법칙이 아니다.

접촉 후보 `zeroLengthContactASDimplex`는 element 자체에 Kn·Kt·μ를 넣는다. μ=tan(접촉마찰각)으로 계산한 값은 JSON에 있다. 현재 시편 구현은 부착력을 포함하지 않으므로 이 논문의 접촉 강도 전체를 재현하지 못한다. 부착력 항이 있는 접촉 모델을 별도로 검증하거나, 무부착 가정을 명시한 비교로 제한해야 한다.

이 실험은 지구에서 수행한 시뮬런트 시험이다. 화성 저중력·온도·기압 전체를 검증한 자료가 아니다. 문헌의 강도 포락선에서 가상 점을 생성해 독립 실험으로 보정하지 않았다.

2026-10-19 원문 재검토: 지반 직접전단은 25/50/75 kPa, 강판 접촉시험은 25/50/75/100 kPa로 구분했다. 기존 normal_stresses_kPa 필드는 지반 직접전단 조건만을 뜻한다. 전단속도·시편 높이·강판 표면 조건을 추가했다. 거친 표면은 에폭시와 모래를 붙인 시험판이므로 일반 강판에 그대로 적용하지 않는다.

건조 시료 압축시험의 하중단계와 문헌 Cc=0.131, Cs=0.030을 추가 기록했다. 이는 체적탄성계수 K나 전단탄성계수 G의 보정값이 아니며, Figure 5 곡선 검토와 별도 압축 시편 검증이 필요하다. 시간 의존 압밀이나 반복진동 검증을 의미하지 않는다.
''',encoding='utf8')
print(json.dumps(rows,indent=2))
