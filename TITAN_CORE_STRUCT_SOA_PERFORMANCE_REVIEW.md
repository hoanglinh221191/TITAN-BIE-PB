# TITAN Core Struct SoA / Performance / Stability Review

Scope: review cac struct va duong du lieu loi quanh `MOL`, `Multi_MOL`, COM, coordinate storage, va muc do tuong thich voi SoA/AVX2/TitanExecutor. Day la review kien truc va rui ro, khong thay doi schema output.

## Findings

### P1 - COM co the chia cho 0 va phat tan NaN/Inf

`MOL::take_com` va `MOL::take_com_all_atoms` tinh COM bang tong khoi luong, nhung chia truc tiep cho `M` khi khong co guard neu tong khoi luong bang 0 hoac residue range rong.

Refs:
- `take_com.cpp:4-48`
- `take_com.cpp:51-76`

Rui ro:
- Neu topology/mass chua duoc gan day du, hoac group selection tao residue khong co atom hop le, COM thanh `nan`/`inf`.
- CMap/QNative mac dinh moi dua tren side-chain COM se bi sai am tham.
- Cac module downstream dung COM de tinh distance/PBC co the sinh matrix/occupancy khong on dinh.

Khuyen nghi:
- Tra ve trang thai thanh cong/that bai cho COM build, hoac ghi `nan` co chu y kem counter/warning.
- Fallback all-atom chi nen chia khi all-atom mass > 0.
- Them smoke cho residue zero-mass, empty side-chain, glycine fallback, va empty range.

### P1 - `Multi_MOL::com_mol` dung topology/parameter khong duoc validate

`Multi_MOL::com_mol` index `parameters[i]` theo so atom cua frame va chia cho `mol_mass`, nhung khong check `num_snap`, `parameters.size()`, hoac `mol_mass > 0`.

Refs:
- `Multi_MOL.cpp:274-285`

Rui ro:
- Frame coordinate va `parameters` lech kich thuoc co the out-of-bounds.
- Tong mass bang 0 tao NaN/Inf.
- Loi co the chi xuat hien voi input/group selection bat thuong, rat kho debug vi fail sau nhieu buoc preprocessing.

Khuyen nghi:
- Them guard `num_snap < m_mol.size()`, `parameters.size() >= m_mol[num_snap].size()`, va `mol_mass > 0`.
- Neu fail, tra error ro rang thay vi tiep tuc tinh.

### P1 - `Multi_MOL::add_coms` co underflow/OOB khi `Index` rong

`add_coms` lap `i < Index.size() - 1` va trong branch con lai dung `Index[0]` ma khong check `Index.empty()`.

Refs:
- `Multi_MOL.cpp:241-263`

Rui ro:
- `Index.size() == 0` gay underflow o dieu kien loop.
- Branch `flag_group_each_atom_as_com` co the doc `Index[0]` out-of-bounds.

Khuyen nghi:
- Guard dau ham: neu `Index.empty()` thi return/error theo semantics hien tai.
- Them test cho group/index rong.

### P2 - Legacy COM distance path copy nguyen `MOL` vao moi thread

`Multi_MOL::cal_dis_res_2_com` tao `std::thread` voi doi so `this->m_mol[snap_num]`, nen moi worker copy ca object `MOL`, gom cac vector toa do lon.

Refs:
- `Multi_MOL.cpp:536-575`

Rui ro:
- Heap churn va memory bandwidth rat lon voi trajectory nhieu frame/atom.
- Copy metadata/string/coordinate lam performance kem va tang nguy co spike RAM.
- Duong nay chua dung `TitanExecutor`, nen khong nhat quan voi huong toi uu Slice D.

Khuyen nghi:
- Toi thieu dung `std::cref(this->m_mol[snap_num])` neu giu raw thread.
- Tot hon: chuyen sang `TitanExecutor::parallel_for`, worker ghi vung rieng, reduce deterministic.
- Validate `Arguments.nt >= 1` hoac normalize thread count.

### P2 - Group extraction co the lam topology va coordinate lech nhau

`take_coordinate_from_group` copy topology theo `base_mol.info[idx]`, sau do voi tung frame lai silently skip coordinate neu `idx >= src.pos.x.size()`.

Refs:
- `take_coordinates.cpp:35-86`

Rui ro:
- `info.size()` cua group frame co the lon hon coordinate count.
- Cac analyzer gia dinh topology/coordinate aligned se doc sai residue/atom hoac out-of-bounds ve sau.

Khuyen nghi:
- Validate group index hop le tren moi frame truoc khi tao grouped molecule.
- Neu muon skip atom, skip ca topology lan coordinate theo mot policy duy nhat va log canh bao.
- Reserve coordinate theo group size de giam realloc.

### P2 - SoA invariant cua `MOL::pos` hien la quy uoc, chua duoc enforce

`MOL::pos` da la SoA (`x`, `y`, `z` vectors rieng), day la nen tot cho SIMD. Tuy nhien cac vector public va co the bi append/clear rieng le tu nhieu call site.

Refs:
- `MOL.h:38-105`
- `MOL.h:203-212`
- `Multi_MOL.cpp:140-147`

Rui ro:
- AVX2/scalar kernels gia dinh `x/y/z` cung size.
- Neu mot call site ghi truc tiep vao `pos.x` ma quen `pos.y/z`, loi co the chay xa moi lo.

Khuyen nghi:
- Them helper invariant `MOL::coordinates_consistent()` va assert/check o bien preprocessing/analyzer.
- Han che ghi truc tiep vao `pos.x/y/z`; uu tien `append_xyz`/builder API.
- Cac kernel SoA nen assert input sizes bang nhau trong debug.

### P2 - `clear_MOL` shrink capacity qua manh cho streaming workload

`MOL::clear_MOL` goi `shrink_to_fit()` cho coordinate va force buffers moi lan clear.

Refs:
- `MOL.h:170-179`

Rui ro:
- Streaming/chunk reuse se mat loi ich reserve, gay allocator churn.
- Khac voi `Multi_MOL::begin_frame`, noi dang clear va reuse capacity tot hon.

Khuyen nghi:
- Tach `clear()` va `release_memory()`.
- Dung shrink chi khi that su can giai phong memory sau run/chunk lon.

### P3 - Force accessors dua vao `assert` nen release build co the UB

`initialize_forces` va force accessors dung `assert(force_vec.has_value())`. Neu build release tat assert, misuse co the thanh undefined behavior.

Refs:
- `MOL.h:149-166`
- `MOL.h:214-239`

Khuyen nghi:
- Doi public accessor sang runtime guard nhe hoac helper `force_enabled()`.
- Giu assert cho debug, nhung release path nen fail ro rang neu caller dung sai.

## SoA Compatibility Assessment

### Diem manh

- `MOL::pos` da dung SoA layout: `x[]`, `y[]`, `z[]`. Day la phan quan trong nhat cho distance kernels, MIC kernels, va AVX2 load lien tiep.
- `get_pos_x/y/z` tra `std::span<const double>`, phu hop voi kernel doc-only va giam copy.
- `titan_contact_kernels` da dat dung huong: tao `PointSoA` rieng cho CMap/QNative, co scalar fallback va AVX2 path.
- `MOL::reduce_forces` da tach hot numeric vectors va co SIMD dispatch cho add-row inplace.

### Diem yeu

- `Multi_MOL` la frame-AoS (`vector<MOL>`). Dieu nay hop voi streaming theo frame, nhung khong phu hop voi vectorization cross-frame. Nen chap nhan cho hien tai vi analyzers dang xu ly frame serial.
- `com` la AoS object co string metadata. Tot cho semantic label/output, nhung khong nen nam tren hot path distance scan.
- CMap/QNative hien phai convert `ContactPoint`/residue point sang temporary SoA moi frame. Dung ve semantics, nhung co allocation/copy overhead.
- Coordinate/topology invariants chua duoc dong goi; SIMD kernels se nhanh nhung can input sach hon de tranh loi am tham.

## Performance Notes

### CMap / QNative

Hien tai huong toi uu Slice D la hop ly:
- build residue/contact points theo semantic cu,
- convert sang SoA,
- scan pair/native distance bang scalar/AVX2 kernel,
- parallel bang `TitanExecutor`.

Nut that con lai:
- per-frame temporary SoA allocation/copy,
- worker-local map merge cua CMap,
- COM build van scalar va phu thuoc metadata/topology.

De toi uu tiep:
- Them scratch buffer tai analyzer level, reuse capacity qua frames.
- Luu `row_index/col_index` va SoA index mapping on dinh khi topology khong doi.
- Tach hot label/key metadata ra khoi hot numeric point arrays.

### Dihedral

Parallel theo residue bang `TitanExecutor` la dung muc cho v1. AVX2 torsion chua dang uu tien vi:
- can `atan2`,
- nhieu branch terminal/chain break/missing atom,
- `nan` handling quan trong hon raw throughput.

### COM

COM build la cau noi giua topology-heavy code va numeric kernels. Neu CMap/QNative chay tren he lon, day co the tro thanh chi phi dang ke.

Huong de nang cap:
- precompute atom index ranges va mass arrays theo residue,
- cache side-chain atom lists/fallback lists,
- tinh COM vao `ResiduePointSoA` truc tiep thay vi tao AoS `com` roi convert lai.

## Stability Checklist Can Them

- `MOL::coordinates_consistent()`: `x.size == y.size == z.size`.
- `MOL::topology_consistent()`: `info.size == pos.size` khi frame co topology.
- `Multi_MOL::parameters_consistent()`: `parameters.size == m_mol[0].info.size`.
- COM guard: zero mass, empty range, invalid residue interval.
- Group extraction guard: group indices hop le tren moi frame.
- `Index` guard: empty/singleton behavior ro rang.

## Implementation Status

Trang thai sau Slice F va pass build/smoke gan nhat:

| Status | Items |
| --- | --- |
| Implemented | COM zero-real-mass fallback dung unit weights/geometric center; `MOL::take_com`, `MOL::take_com_all_atoms`, `Multi_MOL::com_mol`, CMap, va QNative dung cung policy. |
| Implemented | `Multi_MOL::com_mol` guard `num_snap`, coordinate consistency, parameter length, empty frame; `Multi_MOL::add_coms` guard `Index.empty()`, singleton, va invalid `snap_index`. |
| Implemented | `take_coordinate_from_group` validate group indices tren topology va moi frame truoc khi tao grouped molecule; khong con silently tao topology/coordinate lech nhau. |
| Implemented | `MOL::coordinates_consistent()` va `MOL::topology_consistent()` da co; `Multi_MOL::parameters_match_frame()` va `Multi_MOL::analysis_frame_ready()` gom guard topology/coordinate/parameter cho analyzer/core call sites. |
| Implemented | `clear_MOL()` giu capacity cho streaming reuse; `release_memory()` la API ro rang de tra RAM. |
| Implemented | `Multi_MOL::cal_dis_res_2_com` da migrate tu raw `std::thread` sang `TitanExecutor::parallel_for`, moi worker ghi slot rieng va main thread giu `sup_cal_dis_res_com` output path hien tai. |
| Implemented | Force accessors co `MOL::force_enabled()` va runtime guard `std::logic_error("MOL force buffers are not initialized")` khi force buffer chua init. |
| Implemented | CMap/QNative co analyzer-owned scratch buffers cho `ContactPoint` va numeric `PointSoA`, reuse capacity qua frame ma khong doi output schema. |
| Implemented | `titan_contact_geometry` cache dung chung cho CMap/QNative: precompute residue atom lists, side-chain fallback, atom-mode selections, va effective COM weights; per frame chi doc coordinate moi. |
| Implemented | Geometry cache benchmark smoke cho thay cached point build nhanh hon old-style metadata scan tren synthetic 1000 residue fixture trong native smoke run. |
| Future tuning | Neu benchmark he that lon van chi ra COM build la nut that, buoc tiep theo la cache mapping native-reference/trajectory sau hon va scratch-buffer pooling theo worker. |
| Verified by | `compile.bat native all`, `titan_core_struct_smoke`, `titan_contact_geometry_cache_smoke` + benchmark, `titan_contact_kernels_smoke`, `titan_contact_map_smoke`, `titan_qnative_smoke`, `titan_dihedral_smoke`, `cmap_smoke`, `qnative_smoke`, `dihedral_smoke`. |

## Roadmap Completion Status

Core struct/SoA roadmap trong report nay da complete o muc Slice F: cac invariant rui ro cao da co guard, legacy COM distance path da dung `TitanExecutor`, CMap/QNative da dung shared geometry cache va SoA scratch reuse ma khong doi CLI/output/scientific defaults. Phan con lai la future benchmark-driven tuning, khong con la pending blocker cua roadmap nay.

## Suggested Next Slice

1. **Benchmark-driven cache polish**: chi lam neu CMap/QNative tren he lon cho thay geometry build van dang ke sau Slice F.
2. **Worker-local scratch pooling**: can nhac neu profiler thay allocation trong per-frame worker maps/temporary vectors.
3. **Regression discipline**: giu smoke cho zero mass COM, group/index guard, force accessor guard, geometry cache, va CMap/QNative/Dihedral end-to-end.

## Bottom Line

TITAN da co nen SoA kha tot o coordinate layer (`MOL::pos`), CMap/QNative da nam tren duong AVX2/TitanExecutor, va Slice F da dong cac rui ro chinh cua report: COM khong con chia cho 0, topology/coordinate/parameter mismatch co guard ro hon, legacy COM distance path khong con raw-thread copy-heavy, va contact geometry build da duoc cache dung chung. Cac buoc tiep theo nen dua tren profiler he lon thay vi them abstraction truoc.
