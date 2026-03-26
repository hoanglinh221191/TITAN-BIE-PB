#include"Multi_MOL.h"
#include"functions.h"
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<iostream>

int Multi_MOL::topo_protein(string input_topo, string force_f, para& parameters, vector<topo_info_temp> info_topo)
{
	if (input_topo.length() > 5 && input_topo.substr(input_topo.length() - 4) != ".top" && input_topo.substr(input_topo.length() - 4) != ".itp")
	{
		cout << "The file type of topo " << input_topo << " is not supported. This software supports only .top and .itp\n";
		return 1;
	}
	fstream input(input_topo, ios::in);
	if (!input.is_open())
	{
		cout << "Cannot open topo file " << input_topo << "\n";
		return 1;
	}
	else
	{
		cout << "Parsing topo from file " << input_topo << "\n";
	}

	vector<string> content_input;
	
	for (string temp{}; getline(input, temp);)
	{
		if (temp.find(";") != string::npos)
			temp = temp.substr(0,temp.find(";"));
		if (temp.length() < 1)
			continue;
		if (temp.find("#include") != string::npos)
		{
			size_t pos_1{}, pos_2{};
			string temp_2{};
			pos_1 = temp.find("\"");
			if (pos_1 != string::npos)
				temp_2 = temp.substr(pos_1 + 1);
			pos_2=temp_2.find("\"");
			if (pos_2 != string::npos)
			{
				temp_2 = temp.substr(pos_1 + 1, pos_2);
				temp_2.erase(remove(temp_2.begin(), temp_2.end(), ' '), temp_2.end());
				content_input.emplace_back(temp_2);
			}
		}
	}

	if (content_input.empty())
	{
		cout << "cannot parse " << input_topo << " file!!!\n";
		return 1;
	}

	bool flag_in{ false };
	vector<topo_info_temp> info_topo_total_temp;

	for (auto& filename : content_input)
	{
		//cout << filename << "\n";
		string content{}, molname{};
		fstream input(filename, ios::in);
		vector<string> data;
		for (string temp{}; getline(input, temp);)
		{
			/*if (temp.find(";") != string::npos)
				temp = temp.substr(0, temp.find(";"));
			if (temp.length() < 1)
				continue;*/
			data.emplace_back(temp);

			if (temp.find("[ moleculetype ]") != string::npos)
			{
				content = "moleculetype";
			}

			if (content == "moleculetype" && temp.find("[ atoms ]") != string::npos)
			{
				for (size_t j{ data.size() - 1 }; j >= 0; --j)
				{
					if (data[j].find("; name	nrexcl") != string::npos && j+1<data.size())
					{
						molname = temp.substr(0, data[j + 1].find(" "));
						break;
					}
				}
			}

			if (temp.find("[ atomtypes ]") != string::npos)
			{
				content = "atomtypes";
			}
		}

		if (!molname.empty() || input.is_open())
		{
			input.close();
		}

		input.open(filename, ios::in);

		
		for (string temp{}; getline(input, temp);)
		{
			process_atomtypes(input, temp, parameters, force_f, molname);
		}
		input.close();
		int atomnum{};
		input.open(filename, ios::in);
		molname = "Protein";
		vector<topo_info_temp> info_topo_temp;
		for (string temp{}; getline(input, temp);)
		{
			this->process_atom_from_topo(input, temp, molname, atomnum, info_topo_temp);
		}
		if (info_topo_temp.size() == info.size())
		{
			info_topo = info_topo_temp;
			flag_in = true;
		}
		else
		{
			for (auto& i : info_topo_temp)
				info_topo_total_temp.emplace_back(i);
		}
	}
	if (flag_in == false)info_topo = info_topo_total_temp;
	this->post_process_atom_from_topo(info_topo);
	//cout << "info topo size: " << info_topo.size() << '\n';
	//cout << "parameter size: " << parameters.base.size() << '\n';
	return 0;
}