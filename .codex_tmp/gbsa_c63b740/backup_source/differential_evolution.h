#ifndef DE_H
#define DE_H
#include"arguments.h"
#include<sstream>
#include<fstream>
#include<iostream>
#include"miscellaneous.h"

struct de_arg
{
	double diff_weight{ 0.9 }, cross_prob{ 0.9 }, F{ 0.5 }, scan_radius{ 2. }, min_scan_radius{0.5}, restrain_angle{ 65. };
	force old_vector{};
	size_t pop_num{ 2000 }, ite_num{ 100 }, algorithm{};
	bool flag_old_vec{ false };
	de_arg() = default;
	int parse_input(const arguments& Arguments)
	{
		algorithm = Arguments.de_algorithm;
		pop_num = Arguments.de_pop_num;
		ite_num = Arguments.de_ite;
		diff_weight = Arguments.de_diff_weight;
		cross_prob = Arguments.de_cross_prob;
		F = Arguments.de_F;
		scan_radius = fabs(Arguments.de_scan_radius);
		min_scan_radius = fabs(Arguments.de_min_scan_radius);
		restrain_angle = Arguments.de_restrain_angle;
		fstream de_old_vec("DE_out_vector.dat", ios::in);
		if (de_old_vec.is_open())
		{
			double x{}, y{}, z{};
			for (string temp{}; getline(de_old_vec, temp); )
			{
				istringstream is(temp);
				if (is >> x >> y >> z)
				{
					flag_old_vec = true;
					cout << "Found old vector!\n";
					old_vector.Fx = x; old_vector.Fy = y; old_vector.Fz = z;
				}
			}
		}
		return 0;
	}
};
#endif // !DE_H

