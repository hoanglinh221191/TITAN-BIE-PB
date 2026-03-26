#include<string>
#include<string_view>
#include<fstream>
#include<iostream>
#include<algorithm>
#include<sstream>
#include <charconv>
#include <span>
#include <system_error>
#include"Multi_MOL.h"
#include"functions.h"


double fast_parse_double(std::span<const char> data) {
	// Lấy con trỏ bắt đầu và kết thúc từ span
	const char* begin = data.data();
	const char* end = begin + data.size();

	// Bỏ qua khoảng trắng đầu (nếu có)
	while (begin < end && *begin == ' ') {
		++begin;
	}

	double value;
	// from_chars vẫn cần 2 con trỏ
	auto [ptr, ec] = std::from_chars(begin, end, value);

	if (ec == std::errc()) return value;
	return 0.0;
}

int check_and_settime_new(double& time_of_snap, string protein_1, string_view protein_2, arguments& Arguments)
{
	if (protein_1.length() >= 5 && protein_1.find("TITLE") != string::npos)
	{
		size_t found = protein_1.find("=");
		size_t found_2 = protein_1.find("step"), time_length{ protein_1.length() < found_2 ? protein_1.length() : found_2 };

		if (found != string::npos && protein_1.length() > found + 2)
		{
			time_of_snap = (stod(protein_1.substr(found + 1, time_length - found - 1)));
			return 0;
		}
		else
		{
			protein_1 = protein_2;
			found = protein_1.find("t=");
			found_2 = protein_1.find("step");
			time_length = protein_1.length() < found_2 ? protein_1.length() : found_2;
			if (found != string::npos)
			{
				time_of_snap = (stod(protein_1.substr(found + 2, time_length - found - 2)));
			}
			return 0;
		}
	}
	else
		return 1;
}

int Multi_MOL::helper_parse_buffer(vector<string>& buffer_string, arguments& Arguments, box& Box_temp, bool& flag_set_time, double& time_of_snap_temp, size_t& I,
	bool& flag_set_box, ::info& output_info, string input_name, size_t& count_snap, bool& flag_tune_snap, size_t& limit_snap, bool& flag_end)
{
	if (count_snap == limit_snap)
	{
		string_view ref{};
		for (size_t i{}; i < buffer_string.size(); ++i)
		{
			if (i == buffer_string.size() - 1)
				ref = buffer_string[i];
			else
				ref = buffer_string[i + 1];
			if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time, time_of_snap_temp, ref, I, flag_set_box, output_info, buffer_string[i], input_name))
			{
				flag_end = false;
				buffer_string.clear();
				break;
			}
		}
		buffer_string.clear();
		count_snap = 0;

		if (flag_tune_snap == false && !this->m_mol.empty())
		{
			flag_tune_snap = true;
			if (this->m_mol[0].Molecule.size() > 200000)
				limit_snap = 10;
			else if (this->m_mol[0].Molecule.size() > 100000)
				limit_snap = 25;
			else if (this->m_mol[0].Molecule.size() > 30000)
				limit_snap = 70;
			else
				limit_snap = 100;
		}
		return 0;
	}
	else
		return 0;
}

int Multi_MOL::helper_parse_pdb(arguments& Arguments, box& Box_temp, bool& flag_set_time, double& time_of_snap_temp, string_view& protein_2, size_t& I,
	bool& flag_set_box, ::info& output_info, string& protein_1, string input_name)
{

	if (flag_set_time == false)
	{
		time_of_snap_temp = static_cast<double>(I);
		int check_time{ check_and_settime_new(time_of_snap_temp, protein_1, protein_2, Arguments) };
		if (check_time == 0)
			flag_set_time = true;
	}

	if (flag_set_box == false && protein_1.length() >= 54 && protein_1.substr(0, 6) == "CRYST1")
	{
		Box_temp = box(stod(protein_1.substr(6, 9)), stod(protein_1.substr(15, 9)), stod(protein_1.substr(24, 9))
			, stod(protein_1.substr(33, 7)), stod(protein_1.substr(40, 7)), stod(protein_1.substr(47, 7)));
		flag_set_box = true;
	}

	if (protein_1.length() >= 55 && protein_1.substr(0, 4) == "ATOM")
	{
		if (time_of_snap_temp > Arguments.end_time)
			return 1;
		if (Arguments.flag_begin_set == false && I == 0) output_info.set_current_begin_time(time_of_snap_temp);

		if (check_in_time(Arguments, time_of_snap_temp) == true)
		{
			const char* raw = protein_1.data();
			if (I >= 1)
			{
				this->emplace(I, Box_temp, time_of_snap_temp, force_f, mol_info(),
					molecule(fast_parse_double({ raw + 30, 8 }), fast_parse_double({ raw + 38, 8 })
						, fast_parse_double({ raw + 46, 8 })));
			}
			else
			{
				string atom{ protein_1.substr(12, 4) }, resname{ protein_1.substr(16, 5) };
				atom.erase(remove(atom.begin(), atom.end(), ' '), atom.end());
				resname.erase(remove(resname.begin(), resname.end(), ' '), resname.end());

				this->emplace(I, Box_temp, time_of_snap_temp, force_f, mol_info(stoull(protein_1.substr(22, 4)), "NULL", atom, protein_1.substr(21, 1), resname),
					molecule(fast_parse_double({ raw + 30, 8 }), fast_parse_double({ raw + 38, 8 })
						, fast_parse_double({ raw + 46, 8 })));
			}
		}
		else
		{
			//flag_set_time = false;
			flag_set_box = false;
		}
	}
	if ((protein_1 == "ENDMDL" || protein_1 == "END"))
	{
		if (check_in_time(Arguments, time_of_snap_temp) == true)
		{
			if (this->m_mol.empty())
			{
				cout << "Error occured in " << input_name << " file!\nCannot parse atom coordinates of molecule\n";
				return 1;
			}

			check_and_pop_size_2last_snaps(I);

			if (flag_set_box == false && I != 0 && this->m_mol.size() > I)
				this->m_mol[I].Box = this->m_mol[I - 1].Box;
			else
				flag_set_box = false;

			++I;
		}
		flag_set_time = false;
	}
	return 0;
}

int Multi_MOL::parse_pdb(string input_name, string force_f, arguments& Arguments, ::info& output_info)
{
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(NULL);
	fstream input_pdb;
	input_pdb.open(input_name, ios::in);
	if (!input_pdb.is_open())
	{
		cout << "Can not open " << input_name << " file!\n";
		return 1;
	}

	double time_of_snap_temp{};
	box Box_temp;
	bool flag_set_time{ false }, flag_set_box{ false };
	size_t I{}, limit_snap{ 10 };

	bool flag_end{ true }, flag_tune_snap{false};
	while (flag_end)
	{
		vector<string> buffer_string;
		buffer_string.reserve(1000000);
		size_t count_snap{};
		for (string protein_1{}; getline(input_pdb, protein_1);)
		{
			if (protein_1.length() < 1)
				continue;
			buffer_string.emplace_back(move(protein_1));
			if (buffer_string[buffer_string.size()-1] == "ENDMDL" || buffer_string[buffer_string.size() - 1] == "END")
			{
				++count_snap;
				helper_parse_buffer(buffer_string, Arguments, Box_temp, flag_set_time, time_of_snap_temp, I, flag_set_box, output_info, input_name, 
					count_snap, flag_tune_snap, limit_snap, flag_end);
			}
			if (input_pdb.eof())
			{
				flag_end = false;
				break;
			}
			if (flag_end == false)break;
		}

		if (!buffer_string.empty())
		{
			if (buffer_string.back().length() > 1 &&
				buffer_string.back() != "ENDMDL" &&
				buffer_string.back() != "END")
			{
				buffer_string.emplace_back("ENDMDL");
			}
			string_view ref{};
			for (size_t i{}; i < buffer_string.size(); ++i)
			{
				if (i == buffer_string.size() - 1)
					ref = buffer_string[i];
				else
					ref = buffer_string[i + 1];
				if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time, time_of_snap_temp, ref, I, flag_set_box, output_info, buffer_string[i], input_name))
				{
					flag_end = false;
					buffer_string.clear();
					break;
				}
			}
		}
		flag_end = false;
	}

	if (this->m_mol.size() >= 2)
		output_info.set_current_dt(this->m_mol[this->m_mol.size() - 1].time_of_snap - this->m_mol[this->m_mol.size() - 2].time_of_snap);

	if (Arguments.flag_end_set == false && !this->m_mol.empty())
		output_info.set_current_end_time(this->m_mol[this->m_mol.size() - 1].time_of_snap);

	if (Arguments.flag_begin_set == true && this->m_mol.empty() && time_of_snap_temp < Arguments.begin_time)
	{
		cout << "WARNING!!! Your begin time is " << Arguments.begin_time << " ps while time of the last snapshot is smaller than that!\n";
	}
	return 0;
}