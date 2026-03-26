#!/usr/bin/bash
export PATH="/usr/bin:$PATH"
find . -type f \( -name "*.cpp" -o -name "*.h" \) ! -path "*/Eigen/*" ! -path "*/backup*" -exec perl -pi -e 's/push_back\s*\(\s*\{([^}]+)\}\s*\)/emplace_back($1)/g; s/push_back\s*\(\s*(d|od|snap|m|new_mol|val)\s*\)/emplace_back(std::move($1))/g; s/push_back/emplace_back/g;' {} +
