#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>
#include <adios2.h>
#include <numeric>
#include <stdio.h>
#include <math.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
#include <mpi.h>
using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

#if ADIOS2_USE_MPI
    adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD);
#else
    adios2::ADIOS adios("adios2.xml");
#endif

    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " <inputFname> <var L/F> <outputFilename>" << endl;
        MPI_Finalize(); // Ensure MPI is finalized on error
        return -1;
    }
    std::string inputFname = argv[1];
    std::string varIn = argv[2];
    std::string outputFilename = argv[3];

    // making an adios reader
    auto inIO = adios.DeclareIO("WriteIO");
    auto reader = inIO.Open(inputFname, adios2::Mode::Read);
    // creating an engine to write
    adios2::IO bpIO = adios.DeclareIO("time deriv");
    adios2::Engine bpWriter = bpIO.Open(outputFilename, adios2::Mode::Write);

    int step = 0;
    std::vector<double> Fm, Fp, F;

    double deltaT = 0.0;
    int MaxStep = 0;

    std::vector<double> time_deriv;
    std::size_t dataSz = 0;

    adios2::Variable<double> dFdtOut;
    while (true)
    {
        auto status = reader.BeginStep();
        if (status != adios2::StepStatus::OK)
        {
            break;
        }
        int writerStep = reader.CurrentStep();
        std::cout << "Process step " << step
                  << " producer step = " << writerStep << std::endl;

        if (step != writerStep)
        {
            std::cout << "ERROR: this code needs ALL steps from producer to work correctly"
                      << std::endl;
            break;
        }

        auto var = inIO.InquireVariable<double>(varIn);
        // auto shape = var.Shape();

        if (step == 0)
        {
            // reading adios variables
            auto varT = inIO.InquireVariable<double>("deltaT");
            reader.Get<double>(varT, &deltaT, adios2::Mode::Sync);
            std::cout << "   deltaT = " << deltaT << std::endl;
            auto varMaxStep = inIO.InquireVariable<int>("MaxStep");
            reader.Get<int>(varMaxStep, &MaxStep, adios2::Mode::Sync);
        }

        if (step > 0)
        {
            // shifting the values
            std::cout << "    shift data backward" << std::endl;
            Fm = F;
            F = Fp;
        }

        auto shape = var.Shape();
        size_t leny = shape[0];
        size_t lenz = shape[1];
        size_t lenx = shape[2];
        // getting the global size
        dataSz = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<std::size_t>());
        if (step == 0)
        {
            dFdtOut = bpIO.DefineVariable<double>(
                "time_derivative", {leny, lenx, lenz}, {0, 0, 0}, {leny, lenx, lenz}, adios2::ConstantDims);
        }

        Fp.resize(dataSz, 0.0);
        time_deriv.resize(dataSz, 0.0);
        F.resize(dataSz, 0.0);

        reader.Get<double>(var, Fp.data(), adios2::Mode::Sync);

        reader.EndStep();

        if (step >= 2 && step < MaxStep)
        {
            std::cout << "    Calculate with F and DT for " << step - 1 << std::endl;

            // for loop to calculate deriv
            for (std::size_t i = 0; i < dataSz; i++)
            {
                time_deriv[i] = (Fp[i] - Fm[i]) / (2.0 * deltaT);
            }
            bpWriter.BeginStep();
            bpWriter.Put(dFdtOut, time_deriv.data());
            bpWriter.EndStep();
        }
        else if (step == 0)
        {
            time_deriv.resize(dataSz, 0.0);
            bpWriter.BeginStep();
            bpWriter.Put(dFdtOut, time_deriv.data());
            bpWriter.EndStep();
        }

        step++;
    }

    std::fill(time_deriv.begin(), time_deriv.end(), 0.0);
    bpWriter.BeginStep();
    bpWriter.Put(dFdtOut, time_deriv.data());
    bpWriter.EndStep();

    bpWriter.Close();
    MPI_Finalize();

    return 0;
}
