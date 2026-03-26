#ifndef MOL_INFO_H
#define MOL_INFO_H
#include<string>
#include<utility>
#include<cctype>
using namespace std;
class mol_info
{
	size_t resnr, chainnum, respos;
	string atom_type;
	string res_name, original_res_name;
	string atom;
	string chain;
	bool side_chain;
public:
	mol_info() :resnr{  }, atom_type{ "NULL" }, atom{ "NULL" }, chain{ " " }, res_name{ "UND" }, original_res_name{ "UND" }, side_chain{ false }, chainnum{} {}
	template<typename T1, typename T11, typename T12, typename T13>
	mol_info(size_t _resnr, T1&& _atom_type, T11&& _atom, T12&& _chain, T13&& _res_name)
		: resnr{ _resnr }, atom_type{ std::forward<T1>(_atom_type) }, atom{ std::forward<T11>(_atom) }, chain{ std::forward<T12>(_chain) }
		, res_name{ std::forward<T13>(_res_name) }, side_chain{true}
	{
		original_res_name = res_name;
	}
	mol_info(const mol_info& RHS) :resnr{ RHS.resnr }, atom_type{ RHS.atom_type }, atom{ RHS.atom }, chain{ RHS.chain }, res_name{ RHS.res_name }, original_res_name{ RHS.original_res_name }, side_chain{ RHS.side_chain }, chainnum{RHS.chainnum} {}
	mol_info(const mol_info&& RHS) noexcept :resnr{ std::move(RHS.resnr) }, atom_type{ std::move(RHS.atom_type) }, atom{ std::move(RHS.atom) }, chain{ std::move(RHS.chain) }, res_name{ std::move(RHS.res_name) }, 
		original_res_name{ std::move(RHS.original_res_name) }, side_chain{ move(RHS.side_chain) }, chainnum{move(RHS.chainnum)}
		, respos{move(RHS.respos)} {
	}
	mol_info& operator=(const mol_info& rhs)
	{
		resnr = rhs.resnr; atom_type = rhs.atom_type; atom = rhs.atom; chain = rhs.chain;
		res_name = rhs.res_name;
		original_res_name = rhs.original_res_name;
		side_chain = rhs.side_chain;
		chainnum = rhs.chainnum;
		return *this;
	}
	mol_info& operator=(mol_info&& rhs) noexcept
	{
		resnr = std::move(rhs.resnr); atom_type = std::move(rhs.atom_type); atom = std::move(rhs.atom); chain = std::move(rhs.chain);
		res_name = std::move(rhs.res_name);
		original_res_name = std::move(rhs.original_res_name);
		side_chain = move(rhs.side_chain);
		chainnum = move(rhs.chainnum);
		respos = move(rhs.respos);
		return *this;
	}
	~mol_info() = default;
	void set_resnr(size_t res) { resnr = res; }
	const size_t get_resnr() const { return resnr; }
	void set_respos(size_t res) { respos = res; }
	const size_t get_respos()const { return respos; }
	void set_atom_type(const string& type) { atom_type = type; }
	const string get_atom_type() const { return atom_type; }
	void set_chain(const string& cha) { chain = cha; }
	const string get_chain() const { return chain; }
	void set_res_name(const string& res_Name) { res_name = res_Name; }
	void set_res_name(string&& res_Name) { res_name = std::move(res_Name); }
	const string get_res_name() const { return res_name; }
	const string get_original_res_name()const 
	{ 
		if(original_res_name!="HISE") 
			return original_res_name;
	else return "HIS";
	}
	void set_atom(const string& ato) { atom = ato; }
	const string get_atom() const { 
		if (atom.length() >= 4 && isdigit(atom[0]) !=0)
			return atom.substr(1) + atom.substr(0, 1);
		else
			return atom;
	}
	void set_chainnum(size_t num) {
		chainnum = num;
	}
	size_t get_chainnum() const {
		return chainnum;
	}
	bool is_side_chain()const {
		return side_chain;
	}
	bool check_side_chain() {
		if (original_res_name.find("GLY") != string::npos) {
			side_chain = true;
			return true;
		}
		else
			if (atom != "C" && atom != "N" && atom != "CA" && atom != "O"
				&& atom != "HA" && atom != "H" && atom != "H1" && atom != "H2" && atom != "H3")
			{
				side_chain = true;
				return true;
			}
			else
			{
				side_chain = false; return false;
			}
	}
};
#endif
