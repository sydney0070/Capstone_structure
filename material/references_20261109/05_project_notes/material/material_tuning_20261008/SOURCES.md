# 출처와 확인 범위

확인일: 2026-10-08. 논문의 실제 채택 모델과 이번 OpenSeesPy 이식 제안은 [보고서](RESEARCH_REPORT_KO.md)에서 구분했다. 검색 결과의 게시/크롤링 시점은 출판연도로 사용하지 않았다. 원문 전체를 열람하지 못한 자료는 아래에 표시했다.

## 행성 구조물·재료 연구

1. **Mrozek, M., Mrozek, D., Smolana, M., & Anguilano, L. (2025).** *Concept and preliminary structural analysis of a crater-covering dome for future lunar habitats*. Scientific Reports 15, 24744. [출판사 원문](https://www.nature.com/articles/s41598-025-07901-x), [DOI](https://doi.org/10.1038/s41598-025-07901-x). 공개 본문의 Building material, FEM model, boundary conditions 확인. ABAQUS S4R·선형탄성·고정/힌지, SSI 상세 해석 없음. E는 경험식 추정, 휨/직접인장 구분 및 밀도 표기 불일치 기록.
2. **Ruess, F., Schaenzlin, J., & Benaroya, H. (2006).** *Structural Design of a Lunar Habitat*. Journal of Aerospace Engineering 19(3), 133–157. [저자 소속 대학 서지](https://scholarship.libraries.rutgers.edu/esploro/outputs/journalArticle/Structural-Design-of-a-Lunar-Habitat/991031665574004646), [DOI](https://doi.org/10.1061/%28ASCE%290893-1321%282006%2919%3A3%28133%29). 대학 서지·초록과 공개 본문 발췌에서 알루미늄 설계 확인. 전체 지반 코드 재현 근거로 사용하지 않음.
3. **NASA/CR-2005-213845 (2005).** *On Structural Design of a Mobile Lunar Habitat With Multi-Layered Environmental Shielding*. [NASA 원문/서지](https://ntrs.nasa.gov/citations/20050185204). NASA 초록의 탄소복합재 구조 및 폴리에틸렌 차폐층 확인. 보고서의 모든 적층 물성을 코드에 전사하지 않음.
4. **Hybrid lunar inflatable structure (2021).** Acta Astronautica. [출판사](https://www.sciencedirect.com/science/article/abs/pii/S0094576520306342). 공개 초록/발췌의 Aluminum 7075-T6·Kevlar 29 선정 확인. 막 구조는 이번 실행 후보 밖의 확장 대상으로만 기재.
5. **Wan, L., Wendner, R., & Cusatis, G. (2016).** *A Novel Material for In Situ Construction on Mars: Experiments and Numerical Simulations*. Construction and Building Materials 120, 222–231. [저자 공개본](https://arxiv.org/abs/1512.05461), [DOI](https://doi.org/10.1016/j.conbuildmat.2016.05.046). 공개 초록에서 황–시뮬런트, 시험 종류, LDPM 확인. E=15 GPa라는 코드 초기값은 이 논문 추출값이 아님.
6. **Gupta, N., Dawara, V., Kumar, A., & Viswanathan, K. (2023).** *Synthetic space bricks from lunar and martian regolith via sintering*. [저자 공개본](https://arxiv.org/abs/2308.14331). 공개 초록의 LHS/LMS/MGS·소결·압축시험 확인. 코드의 E=20 GPa는 구현용 가정.
7. **On the Structural Design and Additive Construction Process of Martian Habitat Units Using In-Situ Resources on Mars (2025).** Aerospace 12(9), 761. [출판사](https://www.mdpi.com/2226-4310/12/9/761). 검색에 제공된 공개 발췌에서 소재·화성 중력을 적용한 기초 검토 확인. 본문 직접 열람이 제한되어 세부 재료/요소 수치를 실행 코드 근거로 채택하지 않음.

## 지반·시뮬런트

8. **Mechanical Properties of Lunar Regolith and Lunar Soil Simulant.** [NASA 기록 19940004207](https://ntrs.nasa.gov/citations/19940004207). NASA 공개 초록: 입도·밀도·낮은 구속압·시험 경계효과. 기록 번호를 출판연도로 해석하지 않음.
9. **Rame, E., Wilkinson, A., Elliot, A., & Young, C. (2009).** *Flowability of JSC-1a*. NASA/TM-2009-215625. [NASA](https://ntrs.nasa.gov/citations/20090023168). 초록의 전단/다짐 조건별 강도 차이 확인. 모든 JSC-1A 시험 또는 달 현장의 공통 c/φ로 쓰지 않음.
10. **Strength Properties of JSC-1A Lunar Regolith Simulant (2009).** Journal of Geotechnical and Geoenvironmental Engineering 135(5). [DOI/출판사](https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29GT.1943-5606.0000068). 초록의 상대밀도·평균유효응력과 마찰/팽창각 관계 확인; 원시 삼축 곡선은 확보하지 않음.
11. **Development and Geotechnical Engineering Properties of KLS-1 Lunar Simulant (2018).** Journal of Aerospace Engineering 31(1). [출판사](https://ascelibrary.org/doi/10.1061/%28ASCE%29AS.1943-5525.0000798). 공개 초록의 국내 시뮬런트 개발/비교 항목 확인.
12. **Geotechnical Review of Existing Mars Soil Simulants for Surface Mobility.** [NASA 공개 보고서](https://ntrs.nasa.gov/api/citations/20200003046/downloads/20200003046.pdf). 공개 검색 발췌에서 목적별 대표성·역학 자료 문제 확인. 지반 구현의 특정 재료 상수는 전사하지 않음.
13. **Cohesion and shear strength of compacted lunar and Martian regolith simulants (2024).** Icarus. [출판사](https://www.sciencedirect.com/science/article/pii/S0019103524000010), [저자 공개본](https://arxiv.org/abs/2403.11029). 공개 초록과 설명의 밀도·입도·직접전단 관계 확인.
14. **Elgamal, A., Yang, Z., Parra, E., & Ragheb, A. (2003).** *Modeling of cyclic mobility in saturated cohesionless soils*. International Journal of Plasticity 19(6), 883–905. [출판사](https://www.sciencedirect.com/science/article/pii/S0749641902000104), [DOI](https://doi.org/10.1016/S0749-6419%2802%2900010-4). 공개 Model calibration 절에서 Nevada No.120 및 시험 근거 확인. PDMY 계열의 계보로 인용.
15. **Boulanger & Ziotopoulou (2017). PM4Sand v3.1**, UCD/CGM-17/01. [OpenSees 참조/인터페이스](https://openseespydoc.readthedocs.io/en/latest/src/PM4Sand.html), [개발자](https://pm4sand.engr.ucdavis.edu/). OpenSees와 개발자 최신 배포판을 구분.
16. **Boulanger & Ziotopoulou (2018). PM4Silt v1**, UCD/CGM-18/01. [개발자 매뉴얼 목록](https://pm4silt.engr.ucdavis.edu/pm4silt-files), [OpenSees 인터페이스](https://openseespydoc.readthedocs.io/en/latest/src/PM4Silt.html). 공식 설명의 대상·인자·버전 차이 확인.

## SSI·구조 해석 연구와 구현 예제

17. **Boulanger, R. W., Curras, C. J., Kutter, B. L., Wilson, D. W., & Abghari, A. (1999).** *Seismic Soil-Pile-Structure Interaction Experiments and Analyses*. JGGE 125(9), 750–759. [출판사](https://ascelibrary.org/doi/abs/10.1061/%28ASCE%291090-0241%281999%29125%3A9%28750%29). 공개 초록의 점토/모래 원심모형과 동적 p–y 검증 확인.
18. **Harden, C., Hutchinson, T., Martin, G. R., & Kutter, B. L. (2005).** *Numerical Modeling of the Nonlinear Cyclic Response of Shallow Foundations*, PEER 2005-04. [PEER](https://peer.berkeley.edu/publications/2005-04). 공개 초록의 BNWF·모래/점토·모형시험 검증 확인.
19. **Gajan, S., Hutchinson, T. C., Kutter, B. L., Raychowdhury, P., Ugalde, J. A., & Stewart, J. P. (2008).** *Numerical Models for Analysis and Performance-Based Design of Shallow Foundations Subjected to Seismic Loading*, PEER 2007/04, 발행 May 2008. [PEER 보고서](https://peer.berkeley.edu/sites/default/files/web_peer704_sivapalan_gajan_tara_hutchinson_bruce_kutter.pdf). 공개 보고서 검색 텍스트의 기초·스프링 특성 확인. 보고서 번호와 발행연도를 구분.
20. **Raychowdhury & Hutchinson, ShallowFoundationGen documentation (2008).** [공식 PDF](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/documents/ShallowDocumentation.pdf), [공식 명령](https://opensees.berkeley.edu/wiki/index.php/ShallowFoundationGen). BNWF 생성과 수정 스프링 명칭 확인.
21. **McGann & Arduino, Laterally-Loaded Pile Foundation.** [OpenSees 설명](https://opensees.berkeley.edu/wiki/index.php/Laterally-Loaded_Pile_Foundation), [OpenSeesPy 코드](https://openseespydoc.readthedocs.io/en/latest/src/pile.html). 실제 요소·스프링·변환 명령의 구현 참고. 연구 논문과 구분한 예제.
22. **Kramer, S. L., Arduino, P., & Shin, H. (2008).** *Using OpenSees for Performance-Based Evaluation of Bridges on Liquefiable Soils*, PEER 2008-07. [PEER](https://peer.berkeley.edu/publications/2008-07). 공개 초록의 통합 SSI 범위 확인.
23. **Luu, K., Mo, Y. L., & Hsu, T. T. C. (2017).** *Development of CSMM-based shell element for reinforced concrete structures*. Engineering Structures 132, 778–790. [출판사](https://www.sciencedirect.com/science/article/abs/pii/S0141029616314171), [DOI](https://doi.org/10.1016/j.engstruct.2016.11.064). 공개 설명의 CSMM/SCS-3D·OpenSees 기반 구현·검증 대상 확인.
24. **Karakas, A. I., Ozgan, K., & Daloglu, A. T. (2018).** *Soil-structure interaction effects on seismic behavior of a hyperbolic cooling tower using three-parameter Vlasov foundation model*. Earthquakes and Structures 14(1). [출판사 공개 PDF](https://www.techno-press.org/download.php?journal=eas&num=1&ordernum=8&volume=14), [DOI](https://doi.org/10.12989/eas.2018.14.1.085). 실제 해석기는 SAP2000–MATLAB.
25. **Kang, J.-D., Nagae, T., Jeong, S.-H., & Kajiwara, K. (2022).** *Accuracy of Seismic Response Evaluation of Two-Dimensional Analysis Model with Rigid Joints for RC Frame Buildings*. Materials 15(22), 8027. [출판사](https://www.mdpi.com/1996-1944/15/22/8027), [DOI](https://doi.org/10.3390/ma15228027). 공개 검색 발췌에 Concrete02/Steel02 fiber 적용이 나타남. 본문 직접 열람이 제한되어 세부 보정 수치는 채택하지 않음.

## 코드의 공식 명령 근거

26. [PDMY02](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureDependentMultiYield02.html), [PIMY](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/PressureIndependentMultiYield.html), [DruckerPrager](https://opensees.berkeley.edu/wiki/index.php/Drucker_Prager). 명령 인자·재료 단계·출력·DP 압축 자오선 변환.
27. [SSPbrick](https://openseespydoc.readthedocs.io/en/latest/src/SSPbrick.html), [SSPquad](https://openseespydoc.readthedocs.io/en/latest/src/SSPquad.html), [SSPquadUP](https://openseespydoc.readthedocs.io/en/stable/src/SSPquadUP.html). 노드·자유도·출력·동적 u-p 조건.
28. [PM4Sand 공식 단일 요소 개념](https://opensees.berkeley.edu/wiki/index.php/PM4Sand-Example_1), [Python PM4 순환 시험](https://openseespydoc.readthedocs.io/en/latest/src/pm4sand_cyc_cal.html). 경계·단계 전환·압력 DOF 출력 및 Python `setParameter('-val',0,'-ele',1,'FirstCall','1')` 문법 확인. 제공 코드는 자체 재하 경로·출력 체계로 작성. 원 예제의 Rayleigh 감쇠를 복사하지 않음.
29. [Concrete02](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Concrete02.html), [Steel02](https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/Steel02.html). OpenSeesPy Concrete02 페이지의 서두에는 ‘no tensile strength’라는 모순된 문구가 있어 인자 ft/Ets와 OpenSees 본 문서의 선형 인장연화 설명을 따름.
30. [PySimple1](https://openseespydoc.readthedocs.io/en/latest/src/PySimple1.html), [TzSimple1](https://openseespydoc.readthedocs.io/en/latest/src/TzSimple1.html), [QzSimple1](https://openseespydoc.readthedocs.io/en/latest/src/QzSimple1.html). 용량의 힘 단위·type 정의·suction·감쇠.
31. [ElasticPPGap](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/170.htm), [Viscous](https://opensees.berkeley.edu/OpenSees/manuals/usermanual/186.htm), [zeroLengthContact](https://openseespydoc.readthedocs.io/en/latest/src/zeroLengthContact2D.html), [zeroLengthContactASDimplex](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/zeroLengthContactASDimplex.html). 부호·속도 의존성·접촉 인자. 배포 시편에는 마지막 요소의 implicit 적분(`-intType 0`)을 사용하며, 법선 벡터는 2D에서도 3성분을 입력한다.
32. [J2Plasticity](https://openseespydoc.readthedocs.io/en/latest/src/J2Plasticity.html), [ElasticOrthotropic](https://openseespydoc.readthedocs.io/en/latest/src/elasticOrthotropic.html), [PlateFiber wrapper](https://openseespydoc.readthedocs.io/en/latest/src/PlateFiber.html), [PlateFiber section](https://openseespydoc.readthedocs.io/en/latest/src/plateFiberSection.html), [공식 동적 쉘 예제](https://opensees.berkeley.edu/OpenSees/manuals/ExamplesManual/HTML/876.htm).
33. [RC 적층 쉘 공식 Python 예제](https://openseespydoc.readthedocs.io/en/latest/src/RCshearwall.html), [OpenSees 구조 모델링 강습 자료](https://opensees.berkeley.edu/OpenSees/workshops/OpenSeesDays2010/A5_StructModeling.pdf). fiber 보/기둥과 다층 쉘 구현의 구분.
34. [SimpleContact3D](https://opensees.berkeley.edu/wiki/index.php/SimpleContact3D). ContactMaterial3D와의 전용 조합 및 동적 접촉 수치적 주의. 이번 배포에서는 후보 대응표의 확장 경로이며 3D 접촉 시편으로 검증하지 않음.
35. [NASA 행성표](https://nssdc.gsfc.nasa.gov/planetary/factsheet/), [NASA 화성 수치](https://solarsystem.nasa.gov/mars-by-the-numbers/). 중력 근삿값의 근거. 착륙지별 정밀 중력 자료가 아님.

## 자료 처리 원칙

- 블로그·상업 홍보·AI 요약을 물성 보정의 1차 근거로 사용하지 않았다.
- 논문 수치 전체를 복사하지 않았으며 코드의 대부분 초기값과 탐색 범위는 명시적인 예시 가정이다.
- 실험 그래프를 디지타이징한 것처럼 보이는 가상 자료를 만들지 않았다. `synthetic_demo/`는 OpenSees가 생성한 합성 자료임을 기록했다.
- NASA·논문 시뮬런트 시험은 실제 달/화성 현장 시험과 구분했다.
- 공식 문서와 실행 환경의 차이는 실제 실행으로 확인하고 한계를 기록했다.
