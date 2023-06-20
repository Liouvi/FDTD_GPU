#include "dataacc.cuh"
#include <string>
#include <iterator>

void dataacc::extract(std::vector<std::vector<float>> v, std::string name)
{
	std::ofstream output_file(name);
	std::ostream_iterator<float> output_iterator(output_file, "\t");
	for (int i = 0; i < v.size(); i++) {
		copy(v.at(i).begin(), v.at(i).end(), output_iterator);
		output_file << '\n';
	}
}

void dataacc::extract(std::vector<float> v, std::string name)
{

}
