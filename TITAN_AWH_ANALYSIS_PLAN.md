# TITAN AWH Analysis Plan

Tai lieu nay chot roadmap khoa hoc va implementation cho workflow phan tich
GROMACS AWH (Accelerated Weight Histogram) trong Titan. Muc tieu la de lan sau
mo may co the bat dau chay ngay ma khong phai ghep lai y tuong.

================================================================================
1. Summary
================================================================================

Track AWH cua Titan khong nen lay TICA/MSM lam truc chinh. Voi du lieu biased
nhu AWH, workflow trung tam cua v1 la:

- doc reaction coordinate truc tiep tu cac file `awh1.xvg`
- canh hang RC theo thoi gian frame trong trajectory
- chia reaction coordinate thanh cac bin deu nhau
- tim snapshot nam trong tung RC bin
- tinh IFP occupancy theo tung RC bin de mo ta co che unbinding
- xuat representative PDB cho moi bin de lien ket structural interpretation voi RC

Noi cach khac, AWH v1 la `RC-binned structural/state mapping`, khong phai
`TICA-first kinetic model building`.

Nguyen tac dien giai:

- dung ngon ngu `biased ensemble`
- uu tien `structural/state mapping`
- MSM/TPT tren AWH chi duoc xem la `exploratory only`
- khong overclaim `absolute kinetics` hay `equilibrium populations` neu chua co
  reweighting/validation rieng

================================================================================
2. Context and Current Anchors
================================================================================

AWH workflow hien co trong repo va dataset 2336 da co cac diem tu sau:

- multiwalker trajectory list: `awh_walkers.txt`
- single-walker starter: `arg.in`
- multiwalker starter: `arg_multiwalker.in`
- single-walker RC list: `awh_xvg_single.txt`
- multiwalker RC list: `awh_xvg_walkers.txt`
- plotting/report anchor: `plot_titan_awh_multiwalker.py`
- caution/reference text: `discussion_2336_word.html`

Dataset AWH thuc te da xac nhan o:

- `D:\Data_MD\cpps\2336\0\awh1.xtc` ... `D:\Data_MD\cpps\2336\3\awh1.xtc`
- `D:\Data_MD\cpps\2336\0\awh1.xvg` ... `D:\Data_MD\cpps\2336\3\awh1.xvg`

File `awh1.xvg` hien la time series `time(ps) -> position(nm)` cua reaction
coordinate, nen co the duoc dung lam truc phan tich truc tiep ma khong can ep du
lieu biased di qua TICA truoc.

================================================================================
3. Scientific Goals
================================================================================

Track AWH v1 can tra loi duoc cac cau hoi sau:

- ligand roi khoi interface nhu the nao khi di doc biased coordinate?
- nhung contact/H-bond/IFP nao mat truoc, giu lau hon, hoac xuat hien lai theo RC?
- co nhung vung RC nao giong intermediate structural states khong?
- multiwalker co cung phu nhan/ho tro mot co che structural giong nhau khong?

Track AWH v1 khong nham den:

- absolute kinetics
- equilibrium populations duoc dien giai manh
- PMF-native reconstruction
- ket luan TPT/MSM o muc co tinh dinh luong neu chua co reweighting/validation

================================================================================
4. Rollout
================================================================================

4.1 AWH-1 Multiwalker ingest + RC provenance

Status:

- Done trong `titan_awh_ifp.cpp`
- Da co ingest `.xvg` per walker, frame-time matching, matched/unmatched
  accounting, va diagnostics outputs:
  `awh_walker_diagnostics.dat`, `awh_rc_overlap_diagnostics.dat`
- Da verify bang smoke local repo-only duoi:
  `test/awh_inputs_local`, `test/awh_ifp_only_local_native`,
  `test/awh_ifp_local_native`, va `test/awh_ifp_local_portable`

Muc tieu:

- chuan hoa cach doc 1D RC tu `.xvg` cho single-walker va multiwalker
- giu provenance ro rang giua trajectory walker va RC walker

Cong viec:

- doc mot file list `.xvg`, moi dong ung voi mot walker
- parse cot time + RC value tu `.xvg`
- canh hang RC voi frame time trong trajectory da stream
- report so frame match duoc / khong match duoc
- giu log ro rang ve `-b/-e/-dt`, walker count, RC range, va so bin

Deliverables:

- workflow ro rang cho `-traj` va `-i`
- output tom tat RC/bin co provenance day du

4.2 AWH-2 RC binning + representative snapshots

Status:

- Done
- Da co `awh_rc_bin_summary.dat`, `awh_rc_frame_map.dat`, va representative
  PDB theo RC bin
- `awh_rc_frame_map.dat` da duoc enrich them:
  `rc_left`, `rc_center`, `rc_right`, `bin_frames`, `active_walkers`,
  `dominant_walker`, `dominant_fraction`, va `overlap_status`
- Da them `awh_bin_structural_summary.dat` va `awh_bin_top_contacts.csv`
  de doc structural mapping nhanh hon theo tung RC bin

Muc tieu:

- xem AWH ensemble theo coordinate that su duoc bias
- lien ket moi vung RC voi representative structure cu the

Cong viec:

- chia RC thanh cac bin equal-width
- dem so frame trong moi bin
- chon snapshot dai dien gan tam bin nhat
- xuat ban do `bin -> replica_id/time_ps/rc_value`
- xuat representative PDB cho moi bin co du lieu

Deliverables:

- `awh_rc_bin_summary.dat`
- `awh_rc_frame_map.dat`
- bo representative PDB theo RC bin

4.3 AWH-3 IFP along the reaction coordinate

Status:

- Done o muc v1
- Da co `awh_ifp_by_bin.csv` va `awh_ifp_heatmap.dat`
- Da them breakdown tong hop theo loai interaction trong
  `awh_bin_structural_summary.dat`
- Da them top-contact report de rut gon residue-pair noi bat theo bin trong
  `awh_bin_top_contacts.csv`
- Neu can state-transition interpretation giau hon nua, do la phan mo rong
  tiep theo cua AWH-3/AWH-4

Muc tieu:

- dung IFP de mo ta co che unbinding doc theo RC thay vi trong reduced space tru
  tuong

Cong viec:

- gan moi frame vao RC bin
- tinh IFP occupancy theo tung bin
- tong hop contact turnover dọc theo coordinate
- uu tien doc cac tin hieu:
  - contact mat som / mat muon
  - H-bond giu core interface
  - hydrophobic collapse / mo interface
  - residue-pair nao gan voi transition/intermediate RC vung nao

Deliverables:

- `awh_ifp_by_bin.csv`
- `awh_ifp_heatmap.dat`
- cach doc co che unbinding tu IFP occupancy theo RC bin

4.4 AWH-4 Guardrails + optional exploratory analyses

Status:

- Guardrails co ban da done
- Da co warning cho frame khong match, bin trong, va overlap multiwalker kem
- Wording `biased ensemble` / `structural mapping first` da duoc giu o
  roadmap va help
- Da them walker-consistency diagnostics:
  `awh_walker_stability_by_bin.dat` va `awh_walker_pair_overlap.dat`
- Da verify bang fixture consistency local cho thay ro:
  `multiwalker_balanced`, `multiwalker_imbalanced`, va `single_walker`
- Da them provenance `global_frame_index -> RC bin` qua
  `awh_frame_assignments.dat`
- Da them exploratory overlay `awh_exploratory_context.dat` de tong hop
  PCA/UMAP/MSM/TPT context tren biased ensemble khi cac du lieu nay co san
- Overlay nay chi dung cho `exploratory only` state/pathway context;
  khong dien giai unbiased kinetics hay equilibrium populations

Muc tieu:

- giu dien giai khoa hoc an toan tren du lieu biased

Cong viec:

- canh bao khi co qua nhieu frame khong match duoc tu trajectory sang RC `.xvg`
- canh bao khi co bin trong / overlap multiwalker kem
- khuyen nghi so sanh it nhat 2-3 gia tri `-b` de xem ket luan co on dinh khong
- neu dung them PCA/TICA/UMAP/MSM/TPT, coi chung chi la context phu
- gan nhan `exploratory only` cho moi ket luan MSM/TPT tren AWH

Deliverables:

- wording an toan cho report
- checklist guardrail cho biased ensemble interpretation
- `awh_frame_assignments.dat`
- `awh_exploratory_context.dat`

================================================================================
5. Acceptance Criteria
================================================================================

Track nay duoc xem la da du de thuc thi khi:

- single-walker workflow chay ro rang bang `arg.in` + `awh_xvg_single.txt`
- multiwalker workflow chay ro rang bang `arg_multiwalker.in` +
  `awh_xvg_walkers.txt`
- report phat hien duoc truong hop walker overlap kem hoac RC matching kem
- co output IFP theo RC bin, heatmap, frame map, va representative PDBs
- wording MSM/TPT tren AWH mac dinh la `exploratory only`
- co huong dan so sanh `-b` cutoffs de kiem tra do on dinh cua ket luan

================================================================================
6. Assumptions and Guardrails
================================================================================

- Uu tien `structural/state mapping` tren `biased ensemble` truoc.
- Khong xem TICA la truc phan tich mac dinh cho AWH v1.
- Khong xem MSM/TPT la deliverable chinh cho AWH v1.
- Neu can PMF-native diagnostics, histogram-native AWH analysis, hoac
  reweighting-native observables, mo phase rieng ve sau thay vi tron vao v1.
- Khi contact-based co che thay doi manh theo `-b`, uu tien bao cao tinh nhay cam
  nay hon la ep mot ket luan duy nhat.

================================================================================
7. Implementation Starting Points
================================================================================

CLI / workflow surface:

- `-awh-ifp`
- `-awh-xvg <filelist>`
- `-awh-xvg-col <N>`
- `-awh-bins <N>`
- `-awh-out-prefix <prefix>`

Files de bat dau:

- `arg.in`
- `arg_multiwalker.in`
- `awh_walkers.txt`
- `awh_xvg_single.txt`
- `awh_xvg_walkers.txt`
- `plot_titan_awh_multiwalker.py`
- `discussion_2336_word.html`

Output chinh cua implementation hien tai:

- `awh_rc_bin_summary.dat`
- `awh_walker_diagnostics.dat`
- `awh_rc_overlap_diagnostics.dat`
- `awh_walker_stability_by_bin.dat`
- `awh_walker_pair_overlap.dat`
- `awh_ifp_by_bin.csv`
- `awh_ifp_heatmap.dat`
- `awh_bin_structural_summary.dat`
- `awh_bin_top_contacts.csv`
- `awh_rc_frame_map.dat`
- `awh_rep_bin_*.pdb` (hoac prefix tu `-awh-out-prefix`)

Smoke validation inputs da duoc dat duoi:

- `test/awh_ifp_smoke_2336/arg_ifp_only.in`
- `test/awh_ifp_smoke_2336/arg_awh_ifp.in`
- `test/awh_inputs_local/arg_ifp_only_local.in`
- `test/awh_inputs_local/arg_awh_ifp_local.in`
- `test/awh_inputs_local/arg_awh_ifp_consistency_local.in`

Muc dich cua smoke nay la kiem tra:

- workflow AWH RC-binned IFP chay duoc tren dataset 2336
- workflow repo-local AWH smoke chay duoc ma khong phu thuoc `D:\`
- classic `ifp_occupancy.csv` khong doi logic khi bat them `-awh-ifp`
- native / portable giu hash-identical cho toan bo file AWH khoa hoc
- `awh_rc_frame_map.dat` du thong tin hon de doc representative structure theo
  `walker / time / rc / bin`
- co output consistency de chi ro vung RC `single_walker`,
  `multiwalker_balanced`, va `multiwalker_imbalanced`
