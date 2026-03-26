#include"MOL.h"
#include"matrix.h"
#include"groups.h"
#include"functions.h"
#include"folder.h"
#include<fstream>
#include<cmath>
using namespace std;
double gaussian_entropy(const matrix& lamda);
double find_shift(matrix& lamda, int new_row);
void impl_anharmonic(matrix& lamda, matrix& mat_in, arguments& Arguments, int II, int JJ, double& entropy);
double entropy_anharmonic(matrix& lamda, matrix& mat_in, arguments& Arguments);
double entropy_MI_corr(matrix& mat_in, arguments& Arguments);
int matrix::diagonalization(size_t nt)
{
	cout << "Reducing to hessenberg form...\n";
	for (auto i = 0;i < size() - 2;++i)
		reduce_h_form(*this, i + 1, i, nt);
	cout << "Reduce to hessenberg form done\n";
	bool check{ this->check_tol() };
	int ite{1};
	fstream test("test_after_hform.dat", ios::out);
	test << *this;
	int new_row{ row };
	while (!check && ite <= 500 && new_row>1)
	{
		double shift{ find_shift(*this, new_row) };
		if (ite % 100 == 0) {
			cout << ite << "   ";
			cout << "shift " << shift <<"   nRows "<< new_row<< '\n';
		}
		for (int i{};i < row;++i)
			mat[i][i] -= shift;
		new_row=this->given_rotation(new_row, 0,nt);
		for (int i{};i < row;++i)
			mat[i][i] += shift;
		check = this->check_tol();
		++ite;
	}
	cout << "Diagonalization sigma matrix done\n";
	return 0;
}
double imp_cal_entropy(groups_cal& groups, const Multi_MOL& mol1, const Multi_MOL& mol2, arguments& Arguments, size_t group_index)
{
	size_t nt{ Arguments.nt };
	matrix matrix{ static_cast<int>(3 * (groups.groups_atom_number(group_index))), static_cast<int>(mol1.m_mol.size()) };
	for (int j{};j < matrix.col;++j)
		for (int i{};i < (groups.groups_atom_number(group_index));++i)
		{
			size_t atom_index{ groups.atom_index(group_index,i) - 1 };
			if (atom_index>=mol1.m_mol[0].Molecule.size()+ mol2.m_mol[0].Molecule.size())
			{
				cout << "Your entropy group "<< groups.group_name(group_index)<<" is not exist in pdb file!\n";
				exit(EXIT_FAILURE);
			}
			if (atom_index < mol1.m_mol[0].Molecule.size())
			{
				matrix.mat[i * 3][j] = mol1.m_mol[j].Molecule[atom_index].get_x()*sqrt(mol1.parameters[atom_index].get_mass());
				matrix.mat[i * 3 + 1][j] = mol1.m_mol[j].Molecule[atom_index].get_y()*sqrt(mol1.parameters[atom_index].get_mass());
				matrix.mat[i * 3 + 2][j] = mol1.m_mol[j].Molecule[atom_index].get_z()*sqrt(mol1.parameters[atom_index].get_mass());
			}
			else
			{
				size_t mol1size{ mol1.m_mol[0].Molecule.size() };
				matrix.mat[i * 3][j] = mol2.m_mol[j].Molecule[atom_index - mol1size].get_x()*sqrt(mol2.parameters[atom_index - mol1size].get_mass());
				matrix.mat[i * 3 + 1][j] = mol2.m_mol[j].Molecule[atom_index - mol1size].get_y()*sqrt(mol2.parameters[atom_index - mol1size].get_mass());
				matrix.mat[i * 3 + 2][j] = mol2.m_mol[j].Molecule[atom_index - mol1size].get_z()*sqrt(mol2.parameters[atom_index - mol1size].get_mass());
			}

		}
	cout << "Matrix size: "<< matrix.row << " X " << matrix.col << "  \n";
	//cout << "\ngen done\n";
	for (size_t i{};i < matrix.row;++i)
	{
		double average{};
		for (size_t j{};j < matrix.col;++j)
			average += matrix.mat[i][j];
		average /= matrix.col;
		for (size_t j{};j < matrix.col;++j)
			matrix.mat[i][j] -= average;
	}
	::matrix m_temp(matrix);
	::matrix sigma(matrix.row);
	sigma.set_num_threads(Arguments.nt);
	//cout << "gen sigma\n";

	for (size_t i{};i < matrix.col;++i)
	{
		::matrix y_i(matrix.row, 1), y_j(1, matrix.row);
		y_i.set_num_threads(Arguments.nt);
		for (size_t j{};j < matrix.row;++j)
		{
			y_i.mat[j][0] = matrix.mat[j][i];
			y_j.mat[0][j] = matrix.mat[j][i];
		}
		sigma + y_i*y_j;
	}
	//fstream out_mat("out_matrix.dat", ios::out);
	//out_mat << matrix;
	for (size_t i{};i < sigma.row;++i)
		for (size_t j{};j < sigma.col;++j)
			sigma.mat[i][j] /= matrix.col;
	sigma.diagonalization(Arguments.nt);

	//fstream out_sigma("out_sigma.dat", ios::out);
	/*for (size_t i{};i < sigma.row;++i)
	{
		for (size_t j{};j < sigma.col;++j)
		{
			if (fabs(sigma.mat[i][j]) < 0.00005)
				out_sigma << 0.0 << "  ";
			else
				out_sigma << sigma.mat[i][j] << "  ";
		}
		out_sigma << '\n';
	}*/
	double non_gaussian_en{ entropy_anharmonic(sigma,m_temp,Arguments) }, non_gaussian_mi_cor{entropy_MI_corr(m_temp,Arguments)};
	//double non_gaussian_en{ non_gaussian_entropy(m_temp,Arguments.k_neighbor,0) };
	return gaussian_entropy(sigma) - non_gaussian_en - non_gaussian_mi_cor;
}
int cal_entropy(const Multi_MOL& mol1, const Multi_MOL& mol2, groups_cal& groups, arguments& Arguments, 
	 const folders& paths, size_t file_pos, const info& info_output)
{
	if(Arguments.entropy_quasi == true && Arguments.group_entropy_com==-1)
	groups.create_group_entropy(mol1, mol2, Arguments, "calpha");
	double com_entropy{ imp_cal_entropy(groups, mol1, mol2, Arguments, groups.group_entropy_com_index()) };
	cout << "com entropy done, entropy= "<<com_entropy<<"\n----------------------\n";
	double rec_entropy{ imp_cal_entropy(groups, mol1, mol2, Arguments, groups.group_entropy_rec_index()) };
	cout << "rec entropy done, entropy= " <<rec_entropy<<"\n----------------------\n";
	double lg_entropy{ imp_cal_entropy(groups, mol1, mol2, Arguments, groups.group_entropy_lg_index()) };
	cout << "lg entropy done, entropy= "<<lg_entropy<<"\n-----------------------\n";
	string entropy_name{info_output.getEntropyQuasi()};
	if (paths.get_size() > file_pos)
		entropy_name = paths.get_path(file_pos) + "/" + entropy_name;
	fstream entropy_result(entropy_name, ios::out);
	entropy_result << "Entropy (TS) of group " << groups.group_name(groups.group_entropy_com_index()) << ": " << com_entropy << " kJ/mol\n"
		<< "Entropy (TS) of group " << groups.group_name(groups.group_entropy_rec_index()) << ": " << rec_entropy << " kJ/mol\n"
		<< "Entropy (TS) of group " << groups.group_name(groups.group_entropy_lg_index()) << ": " << lg_entropy << " kJ/mol\n"
		<< "-------------------------------------------------------------------\n";
	double max{ com_entropy > rec_entropy ? com_entropy : rec_entropy };
	if (max < lg_entropy) max = lg_entropy;
	entropy_result << "Delta Entropy (T*(S_com-S_rec-S_lg)): " << 2.*max - (com_entropy + rec_entropy + lg_entropy) << " kJ/mol = " << (2.*max - (com_entropy + rec_entropy + lg_entropy)) / 4.184 << " kcal/mol\n";
	entropy_result.close();
	return 0;
}
double gaussian_entropy(const matrix& lamda)
{
	double entropy{}, e{};
	for (size_t i{ 0 };i < lamda.size() -6&& !std::isnan(lamda.mat[i][i]);++i)
	{
		double alpha{ 1.05*10e-13 / sqrt(13.8064 * 300 * fabs(lamda.mat[i][i])) };
		entropy += alpha / (exp(alpha) - 1) - log(1 - exp(-alpha));
		if (std::isinf(log(lamda.mat[i][i])))cout << "[" << i << "]" << "[" << i << "]: " << log(lamda.mat[i][i]) << '\n';
		e += fabs(log(lamda.mat[i][i]));
	}
	return 1.38064*0.3*6.02*entropy;
}
double find_shift(matrix& lamda, int new_row)
{
	double T{ lamda.mat[new_row - 2][new_row - 2] + lamda.mat[new_row - 1][new_row - 1] },
		D{ lamda.mat[new_row - 2][new_row - 2] * lamda.mat[new_row - 1][new_row - 1] - lamda.mat[new_row - 2][new_row - 1] * lamda.mat[new_row - 1][new_row - 2] };
	double E_1{ T / 2. + sqrt(pow(T,2.) / 4.0 - D) }, E_2{ T / 2. - sqrt(pow(T,2.) / 4. - D) };
	if (D < 10e-16 || std::isnan(E_1) || std::isnan(E_2))return lamda.mat[new_row - 1][new_row - 1];
	return(fabs(E_1 - lamda.mat[new_row - 1][new_row - 1]) > fabs(E_2 - lamda.mat[new_row - 1][new_row - 1]) ? E_2 : E_1);
}
void impl_anharmonic(matrix& lamda, matrix& mat_in, arguments& Arguments, int II, int JJ, double& entropy)
{
	for (int i{II};i < JJ && !std::isnan(lamda.mat[i][i]);++i)
	{
		entropy += -1.38064*0.3*6.02*(log(1.05457 / sqrt(13.8064 * 3)) - 0.5* fabs(log(fabs(lamda.mat[i][i]))) + log(sqrt(2 * acos(-1.00)*exp(1))))
			- non_gaussian_entropy(mat_in, Arguments.k_neighbor, i);
	}
}
void imp_mi_corr(matrix& mat_in, arguments& Arguments, int II, int JJ, vector<double>& entropy)
{
	for (int i{ II };i < JJ;++i)
	{
		entropy[i] = non_gaussian_entropy(mat_in, Arguments.k_neighbor, i) / mat_in.row; ///*/ mat_in.row*/+ non_gaussian_entropy(mat_in, Arguments.k_neighbor, j)///* / mat_in.row*/ -

	}
}
double entropy_anharmonic(matrix& lamda, matrix& mat_in, arguments& Arguments)
{
	vector<double> entropy(Arguments.nt,0.0);
	vector<thread> imp(Arguments.nt);
	for (size_t i{};i < Arguments.nt;++i)
		imp[i] = thread(impl_anharmonic, ref(lamda), ref(mat_in), ref(Arguments), static_cast<int>(i*mat_in.row / Arguments.nt), static_cast<int>((i + 1)*mat_in.row / Arguments.nt), ref(entropy[i]));
	for (auto& i : imp)
		i.join();

	double result{};
	for (size_t i{};i < Arguments.nt;++i)
		result += entropy[i];
	cout << "Anharmonic correction done!   "<< result<<"\n";
	return result;
}
void supp_imp_mi_cor(vector<double>& entropy, matrix& mat_in, arguments& Arguments, int I, int J, double& result)
{
	for (int i{ I };i < J;++i)
		for (int j{ i + 1 };j < mat_in.row;++j)
			result += entropy[i] + entropy[j] - non_gaussian_entropy_2nd(mat_in, Arguments.k_neighbor, i, j) / mat_in.row;
}
double entropy_MI_corr(matrix& mat_in, arguments& Arguments)
{
	vector<double> entropy(mat_in.row, 0.0);
	vector<thread> imp(Arguments.nt);
	for (size_t i{};i < Arguments.nt;++i)
		imp[i] = thread(imp_mi_corr, ref(mat_in), ref(Arguments), static_cast<int>(i*mat_in.row / Arguments.nt), static_cast<int>((i + 1)*mat_in.row / Arguments.nt), ref(entropy));
	for (auto& i : imp)
		i.join();

	vector<double> result_ij(Arguments.nt, 0.0);
	for (size_t i{};i < Arguments.nt;++i)
		imp[i] = thread(supp_imp_mi_cor, ref(entropy), ref(mat_in), ref(Arguments), static_cast<int>(i*mat_in.row / Arguments.nt), static_cast<int>((i + 1)*mat_in.row / Arguments.nt), ref(result_ij[i]));
	for (auto& i : imp)
		i.join();
	double result{};
	for (auto& i : result_ij)
		result += i;
	cout << "MI correction done!   "<< result<<"\n";
	return result;
}