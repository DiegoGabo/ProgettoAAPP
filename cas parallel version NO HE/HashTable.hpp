#include <Eigen/Dense>
#include <vector>  
#include <string>

using namespace Eigen;

class HashTable
{
	private:
		int k;
		int m;
		std::vector<unsigned short> thr; //hash+reprobe (16bit)
		std::vector<int> tc; //count (32bit) ..ottimizzazione 2?
		MatrixXf matrix, inverse;
		unsigned short f(unsigned short);
		unsigned short f_rev(unsigned short);
		unsigned short key_vector(std::vector<Nucleotide>);
		int reprobe(unsigned short);
	public:
		HashTable(int); 
		void incrementValue(std::vector<Nucleotide>);
		std::string toString();
};
