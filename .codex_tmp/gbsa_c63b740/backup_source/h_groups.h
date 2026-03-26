#ifndef H_GROUPS_H
#define H_GROUPS_H
#include"MOL.h"
using namespace std;

struct h_group
{
	char name;
	size_t index;
	vector<size_t> h_index;
	h_group() :name{ '-' }, index{} {
		h_index.reserve(100);
	}
	int add_groups(const MOL& mol_in, const MOL& ref, size_t i, const size_t remain, char _name, double bond_length, vector<mol_info>& info);
};

struct h_groups
{
	std::vector<h_group> groups;
	int find_groups(const MOL& mol_in, const MOL& ref, vector<mol_info>& info, const arguments& Arguments);
	void print_group(int m);
};

#endif // !H_GROUPS_H
