#ifndef H_BONDS_H
#define H_BONDS_H
#include<vector>
#include"arguments.h"
using namespace std;
class h_bonds
{
	size_t index_donor, index_acceptor, index_h;
	double angle, dis;
public:
	h_bonds() :index_donor{}, index_acceptor{}, index_h{}, angle{}, dis{} {}
	h_bonds(const h_bonds& rhs)
	{
		this->index_donor = rhs.index_donor;
		this->index_acceptor = rhs.index_acceptor;
		this->index_h = rhs.index_h;
		this->angle = rhs.angle;
		this->dis = rhs.dis;
	}
	template<typename T1, typename T2>
	h_bonds(T1&& _index_donor, T1&& _index_acceptor, T1&& _index_h, T2&& _angle, T2&& _dis) :index_donor{ std::forward<T1>(_index_donor) }
		, index_acceptor{ std::forward<T1>(_index_acceptor) }, index_h{ std::forward<T1>(_index_h) }, angle{ std::forward<T2>(_angle) }, dis{ std::forward<T2>(_dis) } {}
	const size_t get_index_donor() { return index_donor; }
	const size_t get_index_acceptor() { return index_acceptor; }
	const size_t get_index_h()const { return index_h; }
	const double get_angle() { return angle; }
	const double get_dis() { return dis; }
	template<typename T>
	void set_index_donor(T&& i) { index_donor = std::forward<T>(i); }
	template<typename T>
	void set_index_acceptor(T&& i) { index_acceptor = std::forward<T>(i); }
	template<typename T>
	void set_index_h(T&& i) { index_h = std::forward<T>(i); }
	template<typename T>
	void set_angle(T&& i) { angle = std::forward<T>(i); }
	template<typename T>
	void set_dis(T&& i) { dis = std::forward<T>(i); }
	~h_bonds() =default;
};
struct h_bonds_wrapper
{
	vector<h_bonds> dat;
	h_bonds_wrapper() { dat.reserve(200); }
	h_bonds_wrapper(const h_bonds_wrapper& rhs) :dat(rhs.dat) {}
	h_bonds_wrapper& operator=(const h_bonds_wrapper& rhs)
	{
		dat = rhs.dat;
		return *this;
	}
	h_bonds_wrapper& operator=(h_bonds_wrapper&& rhs) noexcept
	{
		dat = move(rhs.dat);
		return *this;
	}
	~h_bonds_wrapper() = default;
};
#endif