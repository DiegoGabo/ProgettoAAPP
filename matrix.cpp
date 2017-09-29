#include <Eigen/Dense>
#include <iostream>
#include <ctime>

using namespace std;
using namespace Eigen;

int main()
{
	int k=3;
	MatrixXf m(2*k, 2*k), inverse(2*k, 2*k);
	srand((unsigned)time(NULL));
	do{
        for(int i=0;i<2*k;i++){
            
			for(int j=0;j<2*k;j++){
                m(i,j) = (rand() % 2);
            	}
        }
	}
	while(m.determinant() == 0);
	std::cout << "\nNormal: " << std::to_string(m.determinant()) << "\n";

	for(int i=0;i<2*k;i++){
         for(int j=0;j<2*k;j++){
            std::cout << std::to_string(m(i,j));
         }
	std::cout << "\n";
	}

	inverse = m.inverse();
	std::cout << "\nInverse: " << std::to_string(inverse.determinant()) << "\n";

	for(int i=0;i<2*k;i++){
         for(int j=0;j<2*k;j++){
            std::cout << std::to_string(inverse(i,j));
         }
	std::cout << "\n";
	}
}
