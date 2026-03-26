#ifndef COM_H
#define COM_H
#include<vector>
#include<string>
#include<utility>
using namespace std;
class com
{
public:
	com() :x{ 0.0 }, y{ 0.0 }, z{ 0.0 }, chain{ " " }, res_name{ "UND" }, res_num{ 0 }, atom_num{}, chainid{} {}
	template<typename T1, typename T2, typename T3, typename T4, typename T5>
	com( T1&& _x, T1&& _y, T1&& _z, T2&& _chain, T3&& _res_name, T4&& _res_num, T4&& _atom_num, T5&& _chainid, T5&& _respos)
		:x{ forward<T1>(_x) }, y{ forward<T1>(_y) }, z{ forward<T1>(_z) }, chain{ forward<T2>(_chain) }, 
		res_name{ forward<T3>(_res_name) }, res_num{ forward<T4>(_res_num) }, atom_num{forward<T4>(_atom_num) }, 
		chainid{ forward<T5>(_chainid) }, respos{forward<T5>(_respos)} {
	}
	~com() = default;
	void set_x(double x_temp){ x = x_temp; }
	double get_x(){ return x; }
	void set_y(double y_temp){ y = y_temp; }
	double get_y(){ return y; }
	void set_z(double z_temp){ z = z_temp; }
	double get_z(){ return z; }
	void set_chain(string chain_temp){ chain = chain_temp; }
	string get_chain(){ return chain; }
	void set_res_name(string res_name_temp){ res_name = res_name_temp; }
	string get_res_name(){ return res_name; }
	void set_res_num(int num){ res_num = num; }
	size_t get_res_num(){ return res_num; }
	void set_atom_num(int _atom_num){ atom_num = _atom_num; }
	size_t get_atom_num(){ return atom_num; }
	void set_chainid(size_t _chainid) { chainid = _chainid; }
	size_t get_chainid() { return chainid; }
	void set_respos(size_t _respos) { respos = _respos; }
	size_t get_respos() { return respos; }
private:
	double x, y, z;
	size_t  res_num, atom_num, chainid, respos;
	string chain,res_name;
};
#endif