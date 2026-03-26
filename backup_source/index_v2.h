#ifndef INDEX_H
#define INDEX_H
#include<string>
#include<vector>
using namespace std;
const size_t N_const = 500;
class index_base
{
	string res_name, type;
	double x, y, z;
	size_t res_pos, res_num;
public:
	index_base() :res_pos{}, res_num{}, x{}, y{}, z{}, res_name{}, type{"unknown"} {}
	index_base(size_t _res_pos) :res_pos{ _res_pos }, x{}, y{}, z{}, res_name{} {}
	template<typename T1, typename T11, typename T2, typename T3>
	index_base(T1&& _res_pos, T11&& _res_num, T2&& _x, T2&& _y, T2&& _z,T3&& _res_name) 
		:res_pos{ forward<T1>(_res_pos) }, res_num{ forward<T11>(_res_num) }, x{ forward<T2>(_x) }, y{ forward<T2>(_y) }, z{ forward<T2>(_z) }, res_name{ forward<T3>(_res_name) } {}
	~index_base() = default;
	index_base(index_base&& rhs) 
		:res_pos{ move(rhs.get_res_pos()) }, res_num{ move(rhs.get_res_num()) }, x{ move(rhs.get_x()) }, y{ move(rhs.get_y()) }, z{ move(rhs.get_z()) }, res_name{ move(rhs.get_res_name()) } {}
	index_base& operator=(const index_base& rhs)
	{
		res_pos = rhs.get_res_pos();
		res_num = rhs.get_res_num();
		x = rhs.get_x();
		y = rhs.get_y();
		z = rhs.get_z();
		res_name = rhs.get_res_name();
		type = rhs.get_type();
		return *this;
	}
	index_base(const index_base& rhs) :res_pos{ rhs.get_res_pos() }, res_num{ rhs.get_res_num() },
		x{ rhs.get_x() }, y{ rhs.get_y() }, z{ rhs.get_z() },
		res_name{ rhs.get_res_name() }, type{rhs.get_type()} {}
	void set_res_pos(size_t _res_pos) { res_pos = _res_pos; }
	void set_res_num(size_t _res_num) { res_num = _res_num; }
	template<typename T>
	inline void set_coor(T&& _x, T&& _y, T&& _z)
	{
		x = forward<T>(_x);
		y = forward<T>(_y);
		z = forward<T>(_z);
	}
	void set_res_name(string& _res_name) { res_name = _res_name; }
	void set_type(string& _type) { type = _type; }
	size_t get_res_pos() const { return res_pos; }
	size_t get_res_num() const { return res_num; }
	double get_x() const { return x; }
	double get_y() const { return y; }
	double get_z() const { return z; }
	string get_res_name() const { return res_name; }
	string get_type() const { return type; }
};
class index
{
	vector<index_base> base;
	string chain{ " " }, type{ "protein" };
	size_t chainid{};
public:
	index() :chain{"NULL"}, type{ "protein" }
	{
		base.reserve(N_const);
	}
	~index() = default;
	template<typename T>
	index(string& _chain, size_t _res_pos, size_t _res_num, T&& _x, T&& _y, T&& _z, string& _res_name) :chain{ _chain }, type{"protein"}
	{
		base.emplace_back( _res_pos, _res_num, _x, _y, _z, _res_name );
	}
	template<typename T>
	index(string& _chain, size_t _chainid, size_t _res_pos, size_t _res_num, T&& _x, T&& _y, T&& _z, string& _res_name) :chain{ _chain }, chainid{ _chainid }, type{ "protein" }
	{
		base.emplace_back(_res_pos, _res_num, _x, _y, _z, _res_name);
	}
	template<typename T1, typename T11, typename T2, typename T3>
	void add_base(T1&& _res_pos, T11&& _res_num, T2&& _x, T2&& _y, T2&& _z, T3&& _res_name)
	{
		base.emplace_back(forward<T1>(_res_pos), forward<T11>(_res_num), forward<T2>(_x), forward<T2>(_y), forward<T2>(_z), forward<T3>(_res_name));
	}
	index(const string& _chain) :chain{ _chain }, type{"protein"} {}
	index(const string& _chain, size_t _chainid) :chain{ _chain }, chainid{ _chainid }, type{ "protein" } {}
	index(string& _chain, size_t _res_pos) :chain{ _chain }, type{ "protein" } { base.emplace_back(index_base(_res_pos)); }
	index(string&& _chain, size_t _res_pos) :chain{ _chain }, type{ "protein" } { base.emplace_back(index_base(_res_pos)); }
	index(string& _chain, size_t _chainid, size_t _res_pos) :chain{ _chain }, chainid{ _chainid }, type{ "protein" } { base.emplace_back(index_base(_res_pos)); }
	index(string&& _chain, size_t _chainid, size_t _res_pos) :chain{ _chain }, chainid{ _chainid }, type{ "protein" } { base.emplace_back(index_base(_res_pos)); }
	size_t get_res_pos_e(size_t _res_pos_e) const { return base[_res_pos_e].get_res_pos(); }
	size_t get_res_num_e(size_t _res_num_e) const { return base[_res_num_e].get_res_num(); }
	string get_chain() const { return chain; }
	size_t get_chainid()const { return chainid; }
	string get_res_name_e(size_t _res_name_e) { return base[_res_name_e].get_res_name(); }
	double get_res_x_e(size_t _res_x_e) const { return base[_res_x_e].get_x(); }
	double get_res_y_e(size_t _res_y_e) const { return base[_res_y_e].get_y(); }
	double get_res_z_e(size_t _res_z_e) const { return base[_res_z_e].get_z(); }
	size_t get_size_v() { return base.size(); }
	void set_res_pos_e(size_t _res_pos_e, size_t _res_pos) { base[_res_pos_e].set_res_pos(_res_pos); }
	void set_res_pos_e(size_t _res_pos) { base.emplace_back(_res_pos, 0u, 0.0, 0.0, 0.0, "NULL"); }
	void set_res_num_e(size_t _res_num_e, size_t _res_num) { base[_res_num_e].set_res_num(_res_num); }
	void set_chain(const string& _chain) { chain = _chain; }
	void set_chainid(const size_t& _chainid) { chainid = _chainid; }
	void set_type(const string& _type) { type = _type; }
	void set_type_base(size_t _base, string type) { base[_base].set_type(type); }
	string get_type() { return type; }
	int get_type_enum() {
		if (type == "dna")
			return 0;
		else if (type == "rna")
			return 1;
		else if (type == "protein")
			return 2;
		else if (type == "complex")
			return 3;
		else if (type == "water")
			return 4;
		else
			return 5;
	}

	string get_type_base(size_t _base) 
	{
	if(_base<base.size()) 
		return base[_base].get_type();
	else 
		return "unknown";
	}

	void set_res_name_e(size_t _res_name_e, string& _res_name) { base[_res_name_e].set_res_name(_res_name); }

	template<typename T>
	inline void set_res_coor_e(const size_t i, T&& _x, T&& _y, T&& _z)
	{
		base[i].set_coor(forward<T>(_x), forward<T>(_y), forward<T>(_z));
	}
};

struct index_molecule
{
	vector<class index> index_mol;
};

#endif
