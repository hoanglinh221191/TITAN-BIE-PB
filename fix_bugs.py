import re

with open("impl_engine.cpp", "r") as f:
    text = f.read()

# Fix double joinable
text = re.sub(r"if \(([a-zA-Z0-9_]+)\.joinable\(\)\)\s*if \(\1\.joinable\(\)\)", r"if (\1.joinable())", text)

# Fix ML string
old_ml = 'if (Arguments.flag_tda)\n          active_in_loop += "TDA ";'
new_ml = 'if (Arguments.flag_tda)\n          active_in_loop += "TDA ";\n        if (Arguments.flag_ml)\n          active_in_loop += "ML ";'
text = text.replace(old_ml, new_ml)

# Move AlignedDoubleVector up
radii_decl = """        AlignedDoubleVector rec_free_radii;
        AlignedDoubleVector lg_free_radii;
        AlignedDoubleVector rec_complex_radii;
        AlignedDoubleVector lg_complex_radii;"""
text = text.replace(radii_decl, "") # remove it from inside

# Add it before the inner loop
text = text.replace("// --- INNER FRAME LOOP ---", 
                    "// --- INNER FRAME LOOP ---\n      AlignedDoubleVector rec_free_radii;\n      AlignedDoubleVector lg_free_radii;\n      AlignedDoubleVector rec_complex_radii;\n      AlignedDoubleVector lg_complex_radii;")

# Replace resize
old_resize = """        if (Arguments.flag_gbsa) {
#ifdef __AVX__
          rec_free_radii.resize(pro_in.m_mol[chunk_snap_idx].size());
          lg_free_radii.resize(active_lg_mol.size());
          rec_complex_radii.resize(pro_in.m_mol[chunk_snap_idx].size());
          lg_complex_radii.resize(active_lg_mol.size());"""

new_resize = """        if (Arguments.flag_gbsa || Arguments.flag_cas) {
          if (rec_free_radii.size() != pro_in.m_mol[chunk_snap_idx].size()) {
            rec_free_radii.resize(pro_in.m_mol[chunk_snap_idx].size());
            rec_complex_radii.resize(pro_in.m_mol[chunk_snap_idx].size());
          }
          std::fill(rec_free_radii.begin(), rec_free_radii.end(), 0.0);
          std::fill(rec_complex_radii.begin(), rec_complex_radii.end(), 0.0);
          
          if (active_lg_mol.size() > 0) {
            if (lg_free_radii.size() != active_lg_mol.size()) {
              lg_free_radii.resize(active_lg_mol.size());
              lg_complex_radii.resize(active_lg_mol.size());
            }
            std::fill(lg_free_radii.begin(), lg_free_radii.end(), 0.0);
            std::fill(lg_complex_radii.begin(), lg_complex_radii.end(), 0.0);
          }
        }

        if (Arguments.flag_gbsa) {
#ifdef __AVX__"""
text = text.replace(old_resize, new_resize)

# Bug 8: Compute radii for CAS if flag_gbsa is false
old_gbsa_end = """          if (t4.joinable())
            t4.join();
#endif"""

new_gbsa_end = """          if (t4.joinable())
            t4.join();
#endif
        } else if (Arguments.flag_cas && !Arguments.flag_gbsa) {
#ifdef __AVX__
          std::thread t1(cal_effective_born_radii_avx2,
                         std::ref(pro_in.m_mol[chunk_snap_idx]),
                         std::ref(rec_vdw_radii), std::ref(rec_free_radii),
                         Arguments.igb, 0, pro_in.m_mol[chunk_snap_idx].size(),
                         nullptr, nullptr);
          std::thread t2;
          std::thread t3;
          std::thread t4;

          if (active_lg_mol.size() > 0) {
            t2 = std::thread(cal_effective_born_radii_avx2,
                             std::ref(active_lg_mol), std::ref(lg_vdw_radii),
                             std::ref(lg_free_radii), Arguments.igb, 0,
                             active_lg_mol.size(), nullptr, nullptr);
            t3 = std::thread(cal_effective_born_radii_avx2,
                             std::ref(pro_in.m_mol[chunk_snap_idx]),
                             std::ref(rec_vdw_radii), std::ref(rec_complex_radii),
                             Arguments.igb, 0, pro_in.m_mol[chunk_snap_idx].size(),
                             &active_lg_mol, &lg_vdw_radii);
            t4 = std::thread(cal_effective_born_radii_avx2,
                             std::ref(active_lg_mol), std::ref(lg_vdw_radii),
                             std::ref(lg_complex_radii), Arguments.igb, 0,
                             active_lg_mol.size(), &pro_in.m_mol[chunk_snap_idx],
                             &rec_vdw_radii);
          }
          if (t1.joinable()) t1.join();
          if (t2.joinable()) t2.join();
          if (t3.joinable()) t3.join();
          if (t4.joinable()) t4.join();
#endif"""
text = text.replace(old_gbsa_end, new_gbsa_end)

with open("impl_engine.cpp", "w") as f:
    f.write(text)

print("Updated impl_engine.cpp")

# --- chunk_streamer.cpp Bug 6 ---
with open("chunk_streamer.cpp", "r") as f:
    chunk_text = f.read()

# Replace assignments at beginning
old_chunk_start = """  pro_in.m_mol.clear();
  lig.m_mol.clear();"""
new_chunk_start = """  if (!pro_in.m_mol.empty()) last_pro_frame = std::move(pro_in.m_mol.back());
  if (!lig.m_mol.empty()) last_lig_frame = std::move(lig.m_mol.back());

  pro_in.m_mol.clear();
  lig.m_mol.clear();"""
chunk_text = chunk_text.replace(old_chunk_start, new_chunk_start)

# Remove the tracking inside eof checks
chunk_text = re.sub(r'      if \(pro_frames > 0\)\n\s*last_pro_frame = pro_in\.m_mol\.back\(\); // Track last known frame\n', '', chunk_text)
chunk_text = re.sub(r'      if \(lig_frames > 0\)\n\s*last_lig_frame = lig\.m_mol\.back\(\); // Track last known frame\n', '', chunk_text)

# Replace the resizing logic at the end
old_resize_logic = """    // Sync shorter trajectory by holding its last known frame
    if (pro_frames < actual_frames) {
      if (pro_frames == 0)
        TitanUI::print_warning(
            "Receptor trajectory ended early. Freezing at last frame.");
      pro_in.m_mol.resize(actual_frames, last_pro_frame);
    }
    if (lig_frames < actual_frames) {
      if (lig_frames == 0)
        TitanUI::print_warning(
            "Ligand trajectory ended early. Freezing at last frame.");
      lig.m_mol.resize(actual_frames, last_lig_frame);
    }"""

new_resize_logic = """    // Sync shorter trajectory by holding its last known frame
    if (pro_frames < actual_frames) {
      if (pro_frames == 0) {
        TitanUI::print_warning(
            "Receptor trajectory ended early. Freezing at last frame.");
        pro_in.m_mol.resize(actual_frames, last_pro_frame);
      } else {
        pro_in.m_mol.resize(actual_frames, pro_in.m_mol.back());
      }
    }
    if (lig_frames < actual_frames) {
      if (lig_frames == 0) {
        TitanUI::print_warning(
            "Ligand trajectory ended early. Freezing at last frame.");
        lig.m_mol.resize(actual_frames, last_lig_frame);
      } else {
        lig.m_mol.resize(actual_frames, lig.m_mol.back());
      }
    }"""

chunk_text = chunk_text.replace(old_resize_logic, new_resize_logic)
with open("chunk_streamer.cpp", "w") as f:
    f.write(chunk_text)

print("Updated chunk_streamer.cpp")
