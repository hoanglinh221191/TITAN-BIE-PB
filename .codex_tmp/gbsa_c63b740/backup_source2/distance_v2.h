#ifndef DISTANCE_H
#define DISTANCE_H
#include<vector>
#include<string>
#include<utility>
using namespace std;
const int N_s = 10000;
struct dmin_base
{
	size_t resnr_rec, resnr_ab, atom_num, atom_num_rec, chainid_rec, chainid_ab, respos_rec;
	double d;
	string chain_rec, chain_ab, res_name_rec, res_name_ab;
	dmin_base() noexcept :resnr_ab{}, resnr_rec{}, atom_num{}, atom_num_rec{}, d{}, chain_ab{}, chain_rec{}, 
		res_name_rec{}, res_name_ab{}, chainid_rec{}, chainid_ab{}, respos_rec{} {
	};

	inline bool operator >(const dmin_base& rhs) const
	{
		return (chainid_rec > rhs.chainid_rec && chainid_ab > rhs.chainid_ab && resnr_rec > rhs.resnr_rec && resnr_ab > rhs.resnr_ab
			&& res_name_rec > rhs.res_name_rec && res_name_ab > rhs.res_name_ab && atom_num > rhs.atom_num
			&& atom_num_rec > rhs.atom_num_rec && respos_rec>rhs.respos_rec);
	}

	inline bool operator < (const dmin_base& rhs) const
	{
		return (chainid_rec < rhs.chainid_rec && chainid_ab < rhs.chainid_ab && resnr_rec < rhs.resnr_rec && resnr_ab < rhs.resnr_ab
			&& res_name_rec < rhs.res_name_rec && res_name_ab < rhs.res_name_ab && atom_num < rhs.atom_num 
			&& atom_num_rec < rhs.atom_num_rec && respos_rec < rhs.respos_rec);
	}

	template<typename _d, typename _chain, typename _chain2, typename _res1, typename _res2, typename _num, typename _num2, typename _num3, typename _num4,
	typename _numchainid, typename _numchainid2>
	dmin_base(_d&& _dis, _chain&& _chain_ab, _chain2&& _chain_rec, _num&& _resnr_ab, _num2&& _resnr_rec,
		_num3&& _atom_num, _num4&& _atom_num_rec, _res1&& _resname_rec, _res2&& _resname_ab, _numchainid&& _chainid_rec, 
		_numchainid2&& _chainid_ab, _numchainid&& _respos_rec)
		:resnr_ab{ _resnr_ab }, resnr_rec{ _resnr_rec }, atom_num{ _atom_num }, atom_num_rec{ _atom_num_rec }, d{ _dis },
		chain_rec{ _chain_rec }, chain_ab{ _chain_ab }, res_name_rec{ _resname_rec }, res_name_ab{ _resname_ab }, chainid_rec{_chainid_rec},
		chainid_ab{ _chainid_ab }, respos_rec{ _respos_rec } {
	}

	inline bool operator == (const dmin_base& rhs) 
	{
		return chainid_rec == rhs.chainid_rec && chainid_ab == rhs.chainid_ab && resnr_ab == rhs.resnr_ab && resnr_rec == rhs.resnr_rec && atom_num == rhs.atom_num && atom_num_rec == rhs.atom_num_rec
			&& res_name_rec == rhs.res_name_rec && res_name_ab == rhs.res_name_ab && respos_rec == rhs.respos_rec;
	}
	~dmin_base() = default;
};
class d_min
{
public:
	d_min() { base.reserve(N_s); }
	vector<dmin_base> base;
	~d_min() = default;
};
#endif
