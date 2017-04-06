#include"cnpy.h"
#include<complex>
#include<cstdlib>
#include<iostream>
#include<map>
#include<string>
#include<vector>

// const int Nx = 3;
// const int Ny = 4;
// const int Nz = 5;

using namespace std;

int main()
{
    vector < vector < float > > genotype_orig(3, vector<float>(4)); // genotype data

    int Nx = (int) genotype_orig.size();
    cout << "n_samples  : " << Nx << "\n";

    int Ny = (int) genotype_orig[0].size();
    cout << "n_genotypes: " << Ny << "\n";

    // generate random data
    cout << "Data:" << "\n";
    float entry;
    for(unsigned i=0;i<Nx;i++){
        for(unsigned j=0;j<Ny;j++){
            entry = rand(); 
            genotype_orig[i][j] = entry;
            cout << entry << "\t";
        }
        cout << "\n";
    }
    
    // reformat data to float vector
    float* data = new float[Nx*Ny];
    for (unsigned i=0;i<Nx;i++){
        for(unsigned j=0;j<Ny;j++){
            data[i*Ny+j] = genotype_orig[i][j];
        }
    }
    const unsigned int shape[] = {Nx,Ny};
    cnpy::npy_save("arr1.npy",data,shape,2,"w");

    // cnpy::NpyArray arr = cnpy::npy_load("arr1.npy");
    // float* loaded_data = reinterpret_cast<float*>(arr.data);

    // for(unsigned k=0;k<Nx*Ny;k++){
    //     cout << data[k] << "\n";
    // }
    delete[] data;
    //create random data
    // std::complex<double>* data = new std::complex<double>[Nx*Ny*Nz];
    // for(int i = 0;i < Nx*Ny*Nz;i++) {
    //     data[i] = std::complex<double>(rand(),rand());
    //     if(i < 10) std::cout << i << ":" << data[i] << '\n'; 
    // }

    // //save it to file
    // const unsigned int shape[] = {Nz,Ny,Nx};
    // cnpy::npy_save("arr1.npy",data,shape,3,"w");

    // //load it into a new array
    // cnpy::NpyArray arr = cnpy::npy_load("arr1.npy");
    // std::complex<double>* loaded_data = reinterpret_cast<std::complex<double>*>(arr.data);
    // 
    // //make sure the loaded data matches the saved data
    // assert(arr.word_size == sizeof(std::complex<double>));
    // assert(arr.shape.size() == 3 && arr.shape[0] == Nz && arr.shape[1] == Ny && arr.shape[2] == Nx);
    // for(int i = 0; i < Nx*Ny*Nz;i++) assert(data[i] == loaded_data[i]);
    // //cleanup: note that we are responsible for deleting all loaded data
    // delete[] data;
    // delete[] loaded_data;
}
