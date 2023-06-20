#ifndef DATAACC_H
#define DATAACC_H
#include <vector>
#include <fstream>
#include <string>


class dataacc
{
public:
	void extract(std::vector<std::vector<float>> v, std::string name);
	void extract(std::vector<float> v, std::string name);
	friend class simsetup;
};

#endif