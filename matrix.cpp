#include <Eigen/Dense>
#include <ctime>

using namespace std;
using namespace Eigen;

int main()
{
	int k = 1; 

	Matrix<int, 2k, 2k> m, inverse;
	

	srand((unsigned)time(NULL));

    do{
        for(int i=0;i<2*k;i++){
            
			for(int j=0;j<2*k;j++){
                m(i,j) = (rand() % 2);
            	}
        }
	}
	while(0);
	cout << "\nNormal: " << "\n";

	for(int i=0;i<2*k;i++){
         for(int j=0;j<2*k;j++){
            cout << std::to_string(m(i,j));
         }
	cout << "\n";
	}
}
    
