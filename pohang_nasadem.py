import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import plotly.graph_objects as go
import webbrowser

import openseespy.opensees as ops


# ============================================================
# 1. 사용자 설정
# ============================================================

# NASADEM HGT 파일 경로
HGT_PATH = Path("n36e129.hgt")

# HGT 타일 정보
# n36e129.hgt는 위도 36~37, 경도 129~130 영역
TILE_LAT_MIN = 36.0
TILE_LON_MIN = 129.0

# 포항 흥해읍 근처 단면 좌표
# 동서 방향 단면
AREA_CENTER_LAT = 36.11
AREA_CENTER_LON = 129.36

# 시각화 영역 크기
AREA_WIDTH_M = 3000.0
AREA_HEIGHT_M = 3000.0
N_X = 50
N_Y = 50

# OpenSees 지반 모델 깊이
BASE_DEPTH = 30.0  # m

# 깊이 방향 요소 개수
N_DEPTH = 12

# 지반 물성값: 1차 단순 모델
E_SOIL = 50.0e6       # Pa = N/m^2
NU_SOIL = 0.30
RHO_SOIL = 1800.0     # kg/m^3

# 모델 두께
THICKNESS = 1.0       # m

# 중력가속도
G = 9.81              # m/s^2


# ============================================================
# 2. NASADEM HGT 읽기 함수
# ============================================================

def read_hgt_file(hgt_path: Path):
    """
    NASADEM/SRTM 계열 .hgt 파일을 읽는다.
    .hgt는 big-endian signed 16-bit integer 형식이다.
    """
    raw = np.fromfile(hgt_path, dtype=">i2")
    n = int(np.sqrt(raw.size))

    if n * n != raw.size:
        raise ValueError(
            f"HGT 파일 크기가 정사각 격자와 맞지 않음: size={raw.size}"
        )

    dem = raw.reshape((n, n)).astype(float)

    # 결측값 처리
    dem[dem == -32768] = np.nan

    return dem


# ============================================================
# 3. 위도/경도 → DEM 인덱스 변환 및 보간
# ============================================================

def latlon_to_pixel(lat, lon, n, lat_min, lon_min):
    """
    HGT 파일은 일반적으로 북쪽에서 남쪽 방향으로 row가 증가한다.
    n36e129의 경우:
      lat: 36~37
      lon: 129~130

    row = 북쪽 37도에서 아래로 내려오는 방향
    col = 서쪽 129도에서 동쪽으로 가는 방향
    """
    lat_max = lat_min + 1.0
    lon_max = lon_min + 1.0

    row = (lat_max - lat) / (lat_max - lat_min) * (n - 1)
    col = (lon - lon_min) / (lon_max - lon_min) * (n - 1)

    return row, col


def bilinear_sample(dem, lat, lon, lat_min, lon_min):
    """
    DEM에서 특정 위도/경도의 고도값을 bilinear interpolation으로 추출한다.
    """
    n = dem.shape[0]
    row, col = latlon_to_pixel(lat, lon, n, lat_min, lon_min)

    r0 = int(np.floor(row))
    c0 = int(np.floor(col))
    r1 = r0 + 1
    c1 = c0 + 1

    if r0 < 0 or c0 < 0 or r1 >= n or c1 >= n:
        return np.nan

    dr = row - r0
    dc = col - c0

    z00 = dem[r0, c0]
    z01 = dem[r0, c1]
    z10 = dem[r1, c0]
    z11 = dem[r1, c1]

    values = np.array([z00, z01, z10, z11], dtype=float)

    if np.any(np.isnan(values)):
        return np.nan

    z = (
        z00 * (1 - dr) * (1 - dc)
        + z01 * (1 - dr) * dc
        + z10 * dr * (1 - dc)
        + z11 * dr * dc
    )

    return z


# ============================================================
# 4. 포항 단면 추출
# ============================================================

def extract_profile(dem):
    """
    포항 흥해읍 근처의 동서 방향 고도 단면을 추출한다.
    """
    lons = np.linspace(PROFILE_LON_START, PROFILE_LON_END, N_PROFILE)
    lats = np.full_like(lons, PROFILE_LAT)

    elevations = np.array([
        bilinear_sample(
            dem,
            lat=lat,
            lon=lon,
            lat_min=TILE_LAT_MIN,
            lon_min=TILE_LON_MIN,
        )
        for lat, lon in zip(lats, lons)
    ])

    # 거리 계산
    # 위도/경도 차이를 근사적으로 m 단위 거리로 변환
    earth_radius = 6371000.0  # m
    lat_rad = np.deg2rad(PROFILE_LAT)

    x = earth_radius * np.cos(lat_rad) * np.deg2rad(lons - lons[0])

    # 결측 제거
    valid = ~np.isnan(elevations)
    x = x[valid]
    lons = lons[valid]
    lats = lats[valid]
    elevations = elevations[valid]

    # 너무 큰 절대 고도보다 상대 고도가 모델링에 편함
    elevations_rel = elevations - np.nanmin(elevations)

    return x, elevations_rel, elevations, lats, lons

# ============================================================
# 4-1. 포항 면적 영역 DEM 추출
# ============================================================

def extract_area_patch(dem):
    """
    중심 위도/경도를 기준으로 일정 크기의 직사각형 DEM 영역을 추출한다.
    결과는 Plotly 3D surface/solid 시각화에 사용한다.
    """

    earth_radius = 6371000.0
    lat_rad = np.deg2rad(AREA_CENTER_LAT)

    # m → degree 변환
    dlat = (AREA_HEIGHT_M / 2.0) / earth_radius * (180.0 / np.pi)
    dlon = (AREA_WIDTH_M / 2.0) / (earth_radius * np.cos(lat_rad)) * (180.0 / np.pi)

    lat_min = AREA_CENTER_LAT - dlat
    lat_max = AREA_CENTER_LAT + dlat
    lon_min = AREA_CENTER_LON - dlon
    lon_max = AREA_CENTER_LON + dlon

    lats = np.linspace(lat_min, lat_max, N_Y)
    lons = np.linspace(lon_min, lon_max, N_X)

    # 실제 거리 좌표계 생성
    x = earth_radius * np.cos(lat_rad) * np.deg2rad(lons - lons[0])
    y = earth_radius * np.deg2rad(lats - lats[0])

    X, Y = np.meshgrid(x, y)

    Z_abs = np.zeros_like(X, dtype=float)

    for iy, lat in enumerate(lats):
        for ix, lon in enumerate(lons):
            Z_abs[iy, ix] = bilinear_sample(
                dem,
                lat=lat,
                lon=lon,
                lat_min=TILE_LAT_MIN,
                lon_min=TILE_LON_MIN,
            )

    # 결측값이 있으면 주변 평균으로 간단 보정
    if np.any(np.isnan(Z_abs)):
        mean_val = np.nanmean(Z_abs)
        Z_abs[np.isnan(Z_abs)] = mean_val

    # 상대 고도
    Z_rel = Z_abs - np.nanmin(Z_abs)

    return X, Y, Z_rel, Z_abs, lats, lons


# ============================================================
# 5. OpenSeesPy 2D 지반 모델 생성
# ============================================================

def build_opensees_soil_model(x, ground_y):
    """
    DEM 지표면을 상부 경계로 하는 2D plane strain 지반 모델 생성.
    """
    ops.wipe()

    # 2D, 각 노드 자유도 ux, uy
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    # 재료 정의
    mat_tag = 1
    ops.nDMaterial("ElasticIsotropic", mat_tag, E_SOIL, NU_SOIL, RHO_SOIL)

    nx = len(x)
    nz = N_DEPTH

    base_y = np.nanmin(ground_y) - BASE_DEPTH

    node_tags = {}

    node_tag = 1

    for i in range(nx):
        for j in range(nz + 1):
            ratio = j / nz

            # j=0: 지표면
            # j=nz: 하부 경계
            y = (1 - ratio) * ground_y[i] + ratio * base_y

            ops.node(node_tag, float(x[i]), float(y))
            node_tags[(i, j)] = node_tag
            node_tag += 1

    # 요소 생성
    ele_tag = 1

    for i in range(nx - 1):
        for j in range(nz):
            # quad 요소 노드 순서: 반시계 방향
            n1 = node_tags[(i, j + 1)]      # bottom-left
            n2 = node_tags[(i + 1, j + 1)]  # bottom-right
            n3 = node_tags[(i + 1, j)]      # top-right
            n4 = node_tags[(i, j)]          # top-left

            ops.element(
                "quad",
                ele_tag,
                n1, n2, n3, n4,
                THICKNESS,
                "PlaneStrain",
                mat_tag,
                0.0,
                RHO_SOIL,
                0.0,
                -G
            )

            ele_tag += 1

    # 경계조건
    # 바닥: x, y 고정
    for i in range(nx):
        bottom_node = node_tags[(i, nz)]
        ops.fix(bottom_node, 1, 1)

    # 좌우 측면: x 방향 고정, y 방향 자유
    for j in range(nz):
        left_node = node_tags[(0, j)]
        right_node = node_tags[(nx - 1, j)]

        ops.fix(left_node, 1, 0)
        ops.fix(right_node, 1, 0)

    return node_tags


# ============================================================
# 6. 고유진동수 해석
# ============================================================

def run_modal_analysis(num_modes=3):
    """
    생성된 OpenSees 모델에 대해 고유치 해석을 수행한다.
    """
    try:
        eigen_values = ops.eigen(num_modes)

        freqs = []

        for lam in eigen_values:
            omega = np.sqrt(lam)
            freq = omega / (2 * np.pi)
            freqs.append(freq)

        return freqs

    except Exception as e:
        print("Modal analysis failed:", e)
        return []


# ============================================================
# 7. 모델 시각화
# ============================================================

def plot_profile(x, ground_y, title="Extracted DEM Profile"):
    plt.figure(figsize=(10, 4))
    plt.plot(x, ground_y, marker="o", markersize=3)
    plt.xlabel("Distance along profile x (m)")
    plt.ylabel("Relative elevation (m)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def plot_mesh(x, ground_y):
    base_y = np.nanmin(ground_y) - BASE_DEPTH
    nz = N_DEPTH

    plt.figure(figsize=(11, 5))

    # 수직선
    for i in range(len(x)):
        ys = []
        xs = []

        for j in range(nz + 1):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            ys.append(y)
            xs.append(x[i])

        plt.plot(xs, ys, linewidth=0.5)

    # 수평선
    for j in range(nz + 1):
        xs = []
        ys = []

        for i in range(len(x)):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            xs.append(x[i])
            ys.append(y)

        plt.plot(xs, ys, linewidth=0.5)

    plt.plot(x, ground_y, linewidth=2, label="DEM ground surface")

    plt.xlabel("Distance x (m)")
    plt.ylabel("Elevation y (m)")
    plt.title("OpenSees 2D Soil Mesh from NASADEM")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
# ============================================================
# 7-4. Plotly HTML 3D 지형 영역 저장
# ============================================================

def save_area_3d_html(X, Y, Z_rel):
    """
    DEM 면적 영역을 3D surface로 시각화한다.
    정사각형 또는 3:4 직사각형 형태의 지형 확인용.
    """

    fig = go.Figure()

    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z_rel,
            colorscale="Earth",
            opacity=0.95,
            colorbar=dict(title="Relative<br>Elevation (m)"),
            name="DEM surface",
        )
    )

    fig.update_layout(
        title="Pohang NASADEM 3D Terrain Area",
        template="plotly_white",
        width=1100,
        height=850,
        scene=dict(
            xaxis_title="East-West distance x (m)",
            yaxis_title="North-South distance y (m)",
            zaxis_title="Relative elevation z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.5, y=-1.8, z=1.1),
            ),
        ),
    )

    html_path = Path("pohang_area_3d_terrain.html").resolve()
    fig.write_html(str(html_path), auto_open=True)

    print(f"Saved: {html_path}")

# ============================================================
# 7-5. Plotly HTML 3D 지반 Solid Block 저장
# ============================================================

def save_area_3d_solid_html(X, Y, Z_rel):
    """
    DEM 면적 영역을 두께가 있는 지반 block처럼 시각화한다.
    Abaqus viewport에서 solid part를 보는 느낌에 가까움.
    """

    base_z = np.nanmin(Z_rel) - BASE_DEPTH

    ny, nx = Z_rel.shape

    x_flat = []
    y_flat = []
    z_flat = []

    top_node = {}
    bottom_node = {}

    node_id = 0

    # top nodes
    for iy in range(ny):
        for ix in range(nx):
            x_flat.append(float(X[iy, ix]))
            y_flat.append(float(Y[iy, ix]))
            z_flat.append(float(Z_rel[iy, ix]))
            top_node[(iy, ix)] = node_id
            node_id += 1

    # bottom nodes
    for iy in range(ny):
        for ix in range(nx):
            x_flat.append(float(X[iy, ix]))
            y_flat.append(float(Y[iy, ix]))
            z_flat.append(float(base_z))
            bottom_node[(iy, ix)] = node_id
            node_id += 1

    I = []
    J = []
    K = []
    intensity = []

    def add_quad(n1, n2, n3, n4, value):
        I.append(n1)
        J.append(n2)
        K.append(n3)
        intensity.append(value)

        I.append(n1)
        J.append(n3)
        K.append(n4)
        intensity.append(value)

    # top surface
    for iy in range(ny - 1):
        for ix in range(nx - 1):
            n1 = top_node[(iy, ix)]
            n2 = top_node[(iy, ix + 1)]
            n3 = top_node[(iy + 1, ix + 1)]
            n4 = top_node[(iy + 1, ix)]

            z_mean = np.mean([z_flat[n1], z_flat[n2], z_flat[n3], z_flat[n4]])
            add_quad(n1, n2, n3, n4, z_mean)

    # bottom surface
    for iy in range(ny - 1):
        for ix in range(nx - 1):
            n1 = bottom_node[(iy, ix)]
            n2 = bottom_node[(iy + 1, ix)]
            n3 = bottom_node[(iy + 1, ix + 1)]
            n4 = bottom_node[(iy, ix + 1)]

            add_quad(n1, n2, n3, n4, base_z)

    # front/back/left/right side surfaces
    # y-min side
    iy = 0
    for ix in range(nx - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = bottom_node[(iy, ix + 1)]
        n3 = top_node[(iy, ix + 1)]
        n4 = top_node[(iy, ix)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n3], z_flat[n4]]))

    # y-max side
    iy = ny - 1
    for ix in range(nx - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = top_node[(iy, ix)]
        n3 = top_node[(iy, ix + 1)]
        n4 = bottom_node[(iy, ix + 1)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n2], z_flat[n3]]))

    # x-min side
    ix = 0
    for iy in range(ny - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = top_node[(iy, ix)]
        n3 = top_node[(iy + 1, ix)]
        n4 = bottom_node[(iy + 1, ix)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n2], z_flat[n3]]))

    # x-max side
    ix = nx - 1
    for iy in range(ny - 1):
        n1 = bottom_node[(iy, ix)]
        n2 = bottom_node[(iy + 1, ix)]
        n3 = top_node[(iy + 1, ix)]
        n4 = top_node[(iy, ix)]
        add_quad(n1, n2, n3, n4, np.mean([z_flat[n3], z_flat[n4]]))

    fig = go.Figure()

    fig.add_trace(
        go.Mesh3d(
            x=x_flat,
            y=y_flat,
            z=z_flat,
            i=I,
            j=J,
            k=K,
            intensity=intensity,
            colorscale="Earth",
            opacity=0.95,
            showscale=True,
            colorbar=dict(title="Elevation<br>(m)"),
            name="Soil solid block",
        )
    )

    # element edge 느낌의 선 추가
    edge_x = []
    edge_y = []
    edge_z = []

    def add_edge(n1, n2):
        edge_x.extend([x_flat[n1], x_flat[n2], None])
        edge_y.extend([y_flat[n1], y_flat[n2], None])
        edge_z.extend([z_flat[n1], z_flat[n2], None])

    # top grid edges
    for iy in range(ny):
        for ix in range(nx - 1):
            add_edge(top_node[(iy, ix)], top_node[(iy, ix + 1)])

    for ix in range(nx):
        for iy in range(ny - 1):
            add_edge(top_node[(iy, ix)], top_node[(iy + 1, ix)])

    # vertical side edges만 일부 표시
    for iy in [0, ny - 1]:
        for ix in range(nx):
            add_edge(top_node[(iy, ix)], bottom_node[(iy, ix)])

    for ix in [0, nx - 1]:
        for iy in range(ny):
            add_edge(top_node[(iy, ix)], bottom_node[(iy, ix)])

    fig.add_trace(
        go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            line=dict(color="black", width=2),
            name="Mesh edges",
            hoverinfo="skip",
        )
    )

    fig.update_layout(
        title="Abaqus-like 3D Solid Terrain Block from NASADEM",
        template="plotly_white",
        width=1200,
        height=850,
        scene=dict(
            xaxis_title="East-West distance x (m)",
            yaxis_title="North-South distance y (m)",
            zaxis_title="Elevation z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.4, y=-1.7, z=1.0),
            ),
        ),
    )

    html_path = Path("pohang_area_3d_solid_abaqus_like.html").resolve()
    fig.write_html(str(html_path), auto_open=True)

    print(f"Saved: {html_path}")
    

# ============================================================
# 8. 단면 좌표 저장
# ============================================================

def save_profile_csv(x, ground_y, elevations_abs, lats, lons):
    data = np.column_stack([x, ground_y, elevations_abs, lats, lons])

    header = "x_m,relative_elevation_m,absolute_elevation_m,latitude,longitude"

    np.savetxt(
        "pohang_nasadem_profile.csv",
        data,
        delimiter=",",
        header=header,
        comments="",
    )

    print("Saved: pohang_nasadem_profile.csv")


# ============================================================
# 9. 메인 실행
# ============================================================

def main():
    print("Reading HGT file...")
    dem = read_hgt_file(HGT_PATH)

    print("DEM shape:", dem.shape)
    print("DEM elevation min:", np.nanmin(dem))
    print("DEM elevation max:", np.nanmax(dem))

    print("Extracting Pohang area patch...")
    X_area, Y_area, Z_area_rel, Z_area_abs, area_lats, area_lons = extract_area_patch(dem)

    print("Area shape:", Z_area_rel.shape)
    print("Area x length:", np.nanmax(X_area) - np.nanmin(X_area), "m")
    print("Area y length:", np.nanmax(Y_area) - np.nanmin(Y_area), "m")
    print("Area elevation range:", np.nanmin(Z_area_rel), np.nanmax(Z_area_rel))

    print("Saving 3D terrain HTML...")
    save_area_3d_html(X_area, Y_area, Z_area_rel)

    print("Saving 3D solid block HTML...")
    save_area_3d_solid_html(X_area, Y_area, Z_area_rel)

    print("Done.")

    #print("Profile points:", len(x))
    #print("Profile length:", x[-1] - x[0], "m")
    #print("Relative elevation range:", np.nanmin(ground_y), np.nanmax(ground_y))

    #save_profile_csv(x, ground_y, elevations_abs, lats, lons)

    # Matplotlib 2D 확인용
    #plot_profile(x, ground_y, title="Pohang Heunghae-eup NASADEM Profile")

    # Plotly HTML 저장
    #save_profile_html(x, ground_y, elevations_abs, lats, lons)

    #print("Building OpenSees model...")
    #build_opensees_soil_model(x, ground_y)

    #print("Number of OpenSees nodes:", len(ops.getNodeTags()))
    #print("Number of OpenSees elements:", len(ops.getEleTags()))

    # Matplotlib 2D mesh 확인용
    #plot_mesh(x, ground_y)

    # Plotly HTML mesh 저장
    #save_mesh_2d_html(x, ground_y)
    #save_mesh_3d_html(x, ground_y)

    #print("Running modal analysis...")
    #freqs = run_modal_analysis(num_modes=3)

    #if freqs:
    #    for i, f in enumerate(freqs, start=1):
    #        print(f"Mode {i}: {f:.4f} Hz")
    #else:
    #    print("Modal analysis result not available.")

    #print("Done.")





# ============================================================
# 10. 시각화
# ============================================================

def save_profile_html(x, ground_y, elevations_abs, lats, lons):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=x,
            y=ground_y,
            mode="lines+markers",
            name="DEM ground surface",
            text=[
                f"lat={lat:.5f}<br>lon={lon:.5f}<br>abs elev={elev:.2f} m"
                for lat, lon, elev in zip(lats, lons, elevations_abs)
            ],
            hovertemplate=(
                "x = %{x:.2f} m<br>"
                "relative elevation = %{y:.2f} m<br>"
                "%{text}<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Pohang Heunghae-eup NASADEM 2D Ground Profile",
        xaxis_title="Distance along profile x (m)",
        yaxis_title="Relative elevation (m)",
        template="plotly_white",
        width=1100,
        height=500,
    )

    html_path = Path("pohang_profile_2d.html").resolve()

    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved: {html_path}")


# ============================================================
# 7-2. Plotly HTML 2D Mesh 저장
# ============================================================

def save_mesh_2d_html(x, ground_y):
    base_y = np.nanmin(ground_y) - BASE_DEPTH
    nz = N_DEPTH

    fig = go.Figure()

    # 수직 mesh line
    for i in range(len(x)):
        xs = []
        ys = []

        for j in range(nz + 1):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            xs.append(x[i])
            ys.append(y)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(width=1, color="gray"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    # 수평 mesh line
    for j in range(nz + 1):
        xs = []
        ys = []

        for i in range(len(x)):
            ratio = j / nz
            y = (1 - ratio) * ground_y[i] + ratio * base_y
            xs.append(x[i])
            ys.append(y)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(width=1, color="gray"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    # 지표면 강조
    fig.add_trace(
        go.Scatter(
            x=x,
            y=ground_y,
            mode="lines+markers",
            line=dict(width=4, color="black"),
            marker=dict(size=4),
            name="DEM ground surface",
        )
    )

    fig.update_layout(
        title="OpenSees 2D Soil Mesh from NASADEM",
        xaxis_title="Distance x (m)",
        yaxis_title="Elevation y (m)",
        template="plotly_white",
        width=1100,
        height=600,
        yaxis=dict(scaleanchor="x", scaleratio=1),
    )

    html_path = Path("pohang_soil_mesh_2d.html").resolve()

    fig.write_html(str(html_path), auto_open=False)

    print(f"Saved: {html_path}")


# ============================================================
# 7-3. Abaqus 느낌의 3D Solid Mesh HTML 저장
# ============================================================

def save_mesh_3d_html(x, ground_y):
    """
    2D plane strain 지반 단면을 시각화용 두께 방향으로 extrusion하여
    Abaqus viewport처럼 보이는 3D solid mesh HTML을 생성한다.

    주의:
    - 해석 모델은 여전히 2D plane strain이다.
    - 여기서 y 방향 두께는 시각화용 가상 두께이다.
    """
    base_z = np.nanmin(ground_y) - BASE_DEPTH
    nx = len(x)
    nz = N_DEPTH

    thickness_vis = 80.0  # 시각화용 두께, 실제 해석 두께 아님
    y_front = -thickness_vis / 2.0
    y_back = thickness_vis / 2.0

    # --------------------------------------------------------
    # 1) 3D 시각화용 노드 생성
    # 좌표계:
    #   X = 단면 진행 방향
    #   Y = 시각화용 두께 방향
    #   Z = 고도/깊이 방향
    # --------------------------------------------------------

    X = []
    Y = []
    Z = []

    front_node = {}
    back_node = {}

    node_id = 0

    for side_name, y_side, node_map in [
        ("front", y_front, front_node),
        ("back", y_back, back_node),
    ]:
        for i in range(nx):
            for j in range(nz + 1):
                ratio = j / nz
                z = (1 - ratio) * ground_y[i] + ratio * base_z

                X.append(float(x[i]))
                Y.append(float(y_side))
                Z.append(float(z))

                node_map[(i, j)] = node_id
                node_id += 1

    # --------------------------------------------------------
    # 2) 표면 face 생성
    # Mesh3d는 삼각형 face만 받으므로 quad를 삼각형 2개로 분할
    # --------------------------------------------------------

    I = []
    J = []
    K = []
    intensity = []

    def add_quad(n1, n2, n3, n4, value):
        # quad: n1-n2-n3-n4
        # triangle 1: n1-n2-n3
        I.append(n1)
        J.append(n2)
        K.append(n3)
        intensity.append(value)

        # triangle 2: n1-n3-n4
        I.append(n1)
        J.append(n3)
        K.append(n4)
        intensity.append(value)

    # front/back face
    for i in range(nx - 1):
        for j in range(nz):
            # front face
            f1 = front_node[(i, j + 1)]
            f2 = front_node[(i + 1, j + 1)]
            f3 = front_node[(i + 1, j)]
            f4 = front_node[(i, j)]

            z_mean_front = np.mean([Z[f1], Z[f2], Z[f3], Z[f4]])
            add_quad(f1, f2, f3, f4, z_mean_front)

            # back face
            b1 = back_node[(i, j + 1)]
            b2 = back_node[(i + 1, j + 1)]
            b3 = back_node[(i + 1, j)]
            b4 = back_node[(i, j)]

            z_mean_back = np.mean([Z[b1], Z[b2], Z[b3], Z[b4]])
            add_quad(b4, b3, b2, b1, z_mean_back)

    # top surface
    for i in range(nx - 1):
        f1 = front_node[(i, 0)]
        f2 = front_node[(i + 1, 0)]
        b2 = back_node[(i + 1, 0)]
        b1 = back_node[(i, 0)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, f2, b2, b1, z_mean)

    # bottom surface
    for i in range(nx - 1):
        f1 = front_node[(i, nz)]
        f2 = front_node[(i + 1, nz)]
        b2 = back_node[(i + 1, nz)]
        b1 = back_node[(i, nz)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, b1, b2, f2, z_mean)

    # left side surface
    for j in range(nz):
        f1 = front_node[(0, j + 1)]
        f2 = front_node[(0, j)]
        b2 = back_node[(0, j)]
        b1 = back_node[(0, j + 1)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, f2, b2, b1, z_mean)

    # right side surface
    for j in range(nz):
        f1 = front_node[(nx - 1, j + 1)]
        f2 = front_node[(nx - 1, j)]
        b2 = back_node[(nx - 1, j)]
        b1 = back_node[(nx - 1, j + 1)]

        z_mean = np.mean([Z[f1], Z[f2], Z[b2], Z[b1]])
        add_quad(f1, b1, b2, f2, z_mean)

    fig = go.Figure()

    # --------------------------------------------------------
    # 3) Solid surface 추가
    # --------------------------------------------------------

    fig.add_trace(
        go.Mesh3d(
            x=X,
            y=Y,
            z=Z,
            i=I,
            j=J,
            k=K,
            intensity=intensity,
            colorscale="Earth",
            opacity=0.92,
            name="Soil solid",
            showscale=True,
            colorbar=dict(title="Elevation<br>(m)"),
            hovertemplate=(
                "x = %{x:.2f} m<br>"
                "y = %{y:.2f} m<br>"
                "z = %{z:.2f} m<br>"
                "<extra></extra>"
            ),
        )
    )

    # --------------------------------------------------------
    # 4) Mesh line 추가
    # Abaqus의 element edge처럼 보이게 하는 부분
    # --------------------------------------------------------

    edge_x = []
    edge_y = []
    edge_z = []

    def add_edge(n1, n2):
        edge_x.extend([X[n1], X[n2], None])
        edge_y.extend([Y[n1], Y[n2], None])
        edge_z.extend([Z[n1], Z[n2], None])

    # front/back mesh lines
    for node_map in [front_node, back_node]:
        # vertical lines
        for i in range(nx):
            for j in range(nz):
                add_edge(node_map[(i, j)], node_map[(i, j + 1)])

        # horizontal lines
        for j in range(nz + 1):
            for i in range(nx - 1):
                add_edge(node_map[(i, j)], node_map[(i + 1, j)])

    # thickness direction connection lines
    for i in range(nx):
        for j in range(nz + 1):
            add_edge(front_node[(i, j)], back_node[(i, j)])

    fig.add_trace(
        go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            line=dict(color="black", width=2),
            name="Element edges",
            hoverinfo="skip",
        )
    )

    # --------------------------------------------------------
    # 5) 지표면 선 강조
    # --------------------------------------------------------

    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=np.full_like(x, y_front),
            z=ground_y,
            mode="lines",
            line=dict(color="red", width=6),
            name="Front ground surface",
        )
    )

    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=np.full_like(x, y_back),
            z=ground_y,
            mode="lines",
            line=dict(color="red", width=6),
            name="Back ground surface",
            showlegend=False,
        )
    )

    # --------------------------------------------------------
    # 6) Layout
    # --------------------------------------------------------

    fig.update_layout(
        title="Abaqus-like 3D Visualization of NASADEM-based OpenSees Soil Mesh",
        template="plotly_white",
        width=1200,
        height=800,
        scene=dict(
            xaxis_title="Distance x (m)",
            yaxis_title="Visualization thickness (m)",
            zaxis_title="Elevation z (m)",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.8, y=-1.8, z=1.1),
                center=dict(x=0.0, y=0.0, z=0.0),
            ),
        ),
        legend=dict(
            x=0.02,
            y=0.98,
        ),
    )

    html_path = Path("pohang_soil_mesh_3d_abaqus_like.html").resolve()

    fig.write_html(str(html_path), auto_open=True)

    print(f"Saved: {html_path}")
    webbrowser.open(html_path.as_uri())
    
if __name__ == "__main__":
    main()