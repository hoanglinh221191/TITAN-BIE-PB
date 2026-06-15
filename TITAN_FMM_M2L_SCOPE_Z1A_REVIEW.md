# Bao Cao Review FMM/M2L Scope Z1A Sau Attribution Z1B

## Ket Luan Da Hieu Chinh

Z1A xac nhan hai ket qua quan trong:

- Tren Barnase, research-medium va research-alt, FMM/M2L bam TABI Direct rat
  sat. Sai so nang luong tuong doi nam trong khoang `2.30e-4` den `7.25e-4`,
  va sampled matvec audit van sach.
- Tren cac he Born/Kirkwood giai tich, FMM va TABI Direct-centroid cho cung ket
  qua nhung cung lech oracle `5-14%` o cac mesh tho. Dieu nay chung minh sai so
  khong den tu ILU/GMRES hay rieng phan far-field FMM.

Nhan dinh ban dau rang centroid near-field la nguyen nhan duy nhat la qua manh.
Z1A chi cho phep xem centroid near quadrature la gia thuyet chinh, vi backend
Direct trong harness cung su dung centroid.

## Bang Chung Z1B

Scope Z1B da chay TABI Direct tren cung geometry voi ba near policies:
`centroid`, `analytic_singular`, va `duffy_adaptive`.

Ma tran gom centered Born va off-center Kirkwood, salt `0.0/0.15 M`, va bon
muc icosphere that `20/80/320/1280` panels. Tat ca 48 rows hoan tat, khong co
near fallback.

Ket qua:

- Tai salt `0.0 M`, ba near modes cho nang luong giong nhau o moi mesh level.
- Tai salt `0.15 M`, non-centroid co ich tren mesh tho, nhung khong group nao
  dat muc cai thien `2x` da dat truoc.
- Cai thien lon nhat la khoang `1.90x` cho off-center Kirkwood 20 panels.
- Tai 1280 panels, cai thien chi con khoang `1.03x` cho Born va `1.12x` cho
  Kirkwood.
- Mesh refinement giam oracle error deu dan cho ca hai fixture. Tai 1280
  panels, error tot nhat con khoang `1.09%` cho Born salt 0.15 M, `0.83%` cho
  Kirkwood salt 0.15 M, va `1.50-1.57%` tai zero salt.

Decision chinh thuc la `z1b_mesh_dominant`.

## Dien Giai Khoa Hoc

Surface discretization/resolution la thanh phan chi phoi cua accuracy gap Z1A.
Near quadrature co dong gop, dac biet voi off-center Yukawa tren mesh tho, nhung
khong phai loi giai duy nhat va loi ich giam khi mesh min hon.

Sampled audit gan machine precision tren cac he nho chi chung minh FMM khop voi
toan tu Direct duoc chon. No khong tu dong chung minh toan tu roi rac khop oracle
vat ly.

## Huong Di Tiep Theo

Khong nen noi FMM-Duffy ngay chi dua tren Z1A. Thu tu hop ly la:

1. Mo rong surface-resolution/discretization gate vuot 1280 panels hoac them
   geometric surface-error metrics.
2. Xac dinh near policy nao con tao loi ich o geometry du chinh xac.
3. Chi sau do moi tao Scope Z2A de noi policy thang vao FMM near path va do
   accuracy/performance tren matched geometry.

Bao cao do luong:
`TITAN_FMM_M2L_SCOPE_Z1B_NEAR_FIELD_ATTRIBUTION_REPORT.txt`.
