#ifndef CACHE_SIGMA_EPSILON_H
#define CACHE_SIGMA_EPSILON_H
#include"miscellaneous.h"
#include"c6c12.h"
#include"arguments.h"
#include<utility>
#include<string>
#include<typeinfo>
#include"force.h"

#ifdef __AVX__
#include <immintrin.h>
#else
#warning AVX is not available.Code will not compile!
#endif

using namespace std;
class cache_force
{
public:
	vector<force> forces;
	cache_force() :forces{} {}
	cache_force(const arguments& Arguments){
		if (Arguments.force_map != 2) forces = vector<force>(2);
	}
	cache_force(cache_force&& rhs): forces(rhs.forces){}
	cache_force(const cache_force& rhs) : forces(rhs.forces) {}

	cache_force(force& _force_elec, force& _force_vdw) {
		forces = vector<force>{ _force_elec,_force_vdw };
	}
	cache_force(force&& _force_elec, force&& _force_vdw) {
		forces = vector<force>{ move(_force_elec), move(_force_vdw) };
	}
	template<typename T>
	cache_force(T&& F_e_x, T&& F_e_y, T&& F_e_z, T&& F_v_x, T&& F_v_y, T&& F_v_z) {
		forces = vector<force>{ force(forward<T>(F_e_x), forward<T>(F_e_y), forward<T>(F_e_z)),
			force(forward<T>(F_v_x) , forward<T>(F_v_y), forward<T>(F_v_z)) };
	}
	cache_force& operator=(const cache_force& rhs);
	cache_force& operator=(cache_force&& rhs);
	cache_force& operator+=(const cache_force& rhs);
	cache_force& operator+=(cache_force&& rhs);

	template<typename T_force>
	void set_force( T_force&& _elec_force, T_force&& _vdw_force)
	{
		if (typeid(T_force) == typeid(force) && !forces.empty())
		{
			forces[0] = forward<T_force>(_elec_force);
			forces[1] = forward<T_force>(_vdw_force);
		}
		else
		{
			forces = vector<force>{ forward<T_force>(_elec_force) ,forward<T_force>(_vdw_force) };
		}
	}
	template<typename T>
	void set_force( T&& elec_fx, T&& elec_fy, T&& elec_fz,
		T&& vdw_fx, T&& vdw_fy, T&& vdw_fz)
	{
		if (!forces.empty())
		{
			forces[0] = force(forward<T>(elec_fx), forward<T>(elec_fy), forward<T>(elec_fz));
			forces[1] = force(forward<T>(vdw_fx), forward<T>(vdw_fy), forward<T>(vdw_fz));
		}
		else
		{
			forces = vector<force>{ force(forward<T>(elec_fx), forward<T>(elec_fy), forward<T>(elec_fz)) ,
				force(forward<T>(vdw_fx), forward<T>(vdw_fy), forward<T>(vdw_fz)) };
		}
	}

	~cache_force() = default;
};
class cache_energy
{
public:
	double elec_ener{}, vdw_ener{};

	cache_energy() : elec_ener{}, vdw_ener{} {}

	cache_energy(const arguments& Arguments) : elec_ener{}, vdw_ener{} {}
	template<typename T>
	cache_energy(T&& _elec, T&& _vdw) : elec_ener{ _elec }, vdw_ener{ _vdw } {}
	template<typename T_elec, typename T_vdw>
	cache_energy(T_elec&& _elec, T_vdw&& _vdw)
		: elec_ener{ _elec }, vdw_ener{ _vdw } {}
	template<typename T_elec, typename T_vdw, typename T>
	cache_energy(T_elec&& _elec, T_vdw&& _vdw): elec_ener{ forward<T>(_elec) }, vdw_ener{ forward<T>(_vdw) } {}

	cache_energy(cache_energy&& rhs) :elec_ener{ move(rhs.elec_ener) }, vdw_ener{ move(rhs.vdw_ener) } {}
	cache_energy(const cache_energy& rhs) :elec_ener{ rhs.elec_ener }, vdw_ener{ rhs.vdw_ener } {}
	cache_energy& operator=(const cache_energy& rhs);
	cache_energy& operator=(cache_energy&& rhs);
	cache_energy& operator+=(const cache_energy& rhs);
	cache_energy& operator+=(cache_energy&& rhs);

	~cache_energy() = default;

	template<typename T>
	void set_elec_vdw(T&& elec, T&& vdw)
	{
		elec_ener = std::forward<T>(elec);
		vdw_ener = std::forward<T>(vdw);
	}
	template<typename T_elec, typename T_vdw>
	void set_energy(T_elec&& _elec_ener, T_vdw&& _vdw_ener)
	{
		elec_ener = forward<T_elec>(_elec_ener);
		vdw_ener = forward<T_vdw>(_vdw_ener);
	}

};

struct cached_value
{
	double sigma, epsilon;
	cached_value() :sigma(), epsilon() {}
	template<typename T1, typename T2>
	cached_value(T1&& _sigma, T2&& _epsilon) : sigma{ forward<T1>(_sigma) }, epsilon{ forward<T2>(_epsilon) } {}
};

struct cache_sigma_epsilon
{
	vector<cache_energy> parameters;
	vector<cache_force> force_cache;
	vector<double> charge;
	string force_field_name;
	double cutoff_dis{ -1.1 };
	cache_sigma_epsilon() = default;
	cache_sigma_epsilon(double dis_cutoff) :cutoff_dis{ dis_cutoff }, force_field_name{"amber99sb"} 
	{
		parameters.reserve(10000000);
		force_cache.reserve(1000000);
	}
	cache_sigma_epsilon(size_t MOL_size) : force_field_name{"amber99sb"}, 
		parameters{ vector<cache_energy>(MOL_size) },force_cache{ vector<cache_force>(MOL_size) }{}

	cache_sigma_epsilon(size_t MOL_size, const arguments& Arguments): force_field_name{"amber99sb"}, parameters{ vector<cache_energy>(MOL_size, Arguments) }
	{
		if (Arguments.flag_cal_force == 1)
			force_cache = vector<cache_force>(MOL_size, Arguments);
	}

	int set_size(size_t MOL_size, const arguments& Arguments, size_t mol_1_size, size_t mol_2_size);
	void set_elec_vdw(size_t pos, double* elec, double* vdw, size_t sz)
	{
		for (size_t i{}; i < sz; ++i)
			parameters[pos + i].set_elec_vdw(elec[i], vdw[i]);
	}

};
#endif // !CACHE_SIGMA_EPSILON_H