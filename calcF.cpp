
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>
#include <adios2.h>
#include <stdio.h>
#include <math.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
#include <mpi.h>
#define PI 3.141592653589793
using namespace std;
// when the input number get larger Laplace is wrong
void gencore(int rank, double h, size_t start_y, std::vector<double> &x, std::vector<double> &z, std::vector<double> &y)
{

    for (int i = 0; i < x.size(); i++)
    {
        x[i] = (i * h);
    }
    for (int k = 0; k < z.size(); k++)
    {
        z[k] = (k * h);
    }

    for (int j = 0; j < y.size(); j++)
    {

        y[j] = ((start_y + j) * h);
    }
}

void evalF(int num, int rank, double t, std::vector<double> &x, std::vector<double> &z, std::vector<double> &y, std::vector<double> &F)
{
    double theta = fmod(t, 2.0 * PI);
    double phi = fmod(t / 2.0, 2.0 * PI); // for nums
    double x0 = sin(phi / 4.0) + 3.0;
    double y0 = sin(phi / 4.0) + 3.0;
    double z0 = PI;
    double spin = -2.0;
    double x1 = spin * sin(phi) + x0;
    double y1 = spin * cos(phi) + y0;
    double z1 = PI;
    double scalar = 2 * PI * 2 * PI * 2 * PI;

    for (int j = 0; j < y.size(); j++) //  j= row
    {
        for (int k = 0; k < z.size(); k++)
        {
            for (int i = 0; i < x.size(); i++) // i =col
            {
                int L = i + x.size() * k + j * x.size() * z.size();
                if (num == 1) // boring
                {
                    F[L] = exp((x[i] * x[i] + y[j] * y[j] + z[k] * z[k]) / scalar);
                }
                else if (num == 2) // cool
                {
                    F[L] = 1.0 / (sqrt(pow(x[i] - x0, 2) + pow(y[j] - y0, 2) + pow(z[k] - z0, 2)));
                }
                else if (num == 3) // intresting
                {
                    F[L] = 1.0 / (sqrt(pow(x[i] - x0, 2) + pow(y[j] - y0, 2) + pow(z[k] - z0, 2))) + 1.0 / (sqrt(pow(x[i] - x1, 2) + pow(y[j] - y1, 2) + pow(z[k] - z1, 2)));
                }
                else // WOWERS
                {
                    F[L] = 1.5 / (pow((x[i] - x0) * cos(theta) - (y[j] - y0) * sin(theta), 2) + 4 * pow((x[i] - x0) * sin(theta) + (y[j] - y0) * cos(theta), 2) + pow(z[k] - z0, 2)) + 0.5 / (pow((x[i] - x1) * cos(theta) - (y[j] - y1) * sin(theta), 2) + 4 * pow((x[i] - x1) * sin(theta) + (y[j] - y1) * cos(theta), 2) + pow(z[k] - z1, 2));
                }
            }
        }
    }
}

int main(int argc, char *argv[])

{
    int rank, size;
#if ADIOS2_USE_MPI
    // std::cout << "Enter file name, y size, z size, x size, time steps, function number 1-4" << endl;
    int provided;
    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    // rank is giving a number to each proccess
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //  the size of something
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#else
    // std::cout << "Enter file name, y size, z size, x size, time steps, function number 1-4" << endl;
    rank = 0;
    size = 1;
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD);
#else
    adios2::ADIOS adios("adios2.xml");
#endif

    if (argc < 7)
    {
        cout << "Minimum 6 arguments required!!!!" << std::endl;
        cout << "ethan <filename> <global_ny> <global_nz> <global_nx> <steps> <which equation: 1=x^2+y^2+z^2" << endl;
        return 1;
    }
    std::string filename;
    filename = argv[1];
    size_t lenx = atoi(argv[2]);
    size_t leny = atoi(argv[3]);
    size_t lenz = atoi(argv[4]);
    size_t nsteps = atoi(argv[5]);
    int num = atoi(argv[6]);

    size_t local_len_y = leny / size;
    size_t remainder = leny % size;
    if (rank < remainder)
    {
        local_len_y++;
    }
    size_t len_global = lenx * leny * lenz;
    size_t len_local = local_len_y * lenx * lenz;

    size_t start_y = rank * (leny / size) + (rank < remainder ? rank : remainder);

    std::vector<double> x, y, z, tarr;
    x.resize(lenx, 0);
    z.resize(lenz, 0);
    y.resize(local_len_y, 0);

    double dt = 1.0 / nsteps;
    double h = 2 * PI / (lenx - 1);
    std::vector<double> F(len_local);
    F.resize(len_local, 0.0);

    std::vector<double> laplace(len_local);
    laplace.resize(len_local, 0.0);

    adios2::IO bpIO = adios.DeclareIO("WriteIO");
    adios2::Engine bpWriter = bpIO.Open(filename, adios2::Mode::Write);
    adios2::Operator op = adios.DefineOperator("mgard", "mgard");

    adios2::Variable<double> xOut = bpIO.DefineVariable<double>(
        "x", {lenx}, {0}, {lenx}, adios2::ConstantDims);
    adios2::Variable<double> zOut = bpIO.DefineVariable<double>(
        "z", {lenz}, {0}, {lenz}, adios2::ConstantDims);
    adios2::Variable<double> yOut = bpIO.DefineVariable<double>(
        "y", {leny}, {start_y}, {local_len_y}, adios2::ConstantDims);

    adios2::Variable<double> fOut = bpIO.DefineVariable<double>(
        "F", {leny, lenz, lenx}, {start_y, 0, 0}, {local_len_y, lenz, lenx}, adios2::ConstantDims);
    // SOMETHING THAT USED TO WORK WITH CMAKE??? BUT WORKS WITH BEING UNCOMMENTED BUT THE XML IS NOT? IDK
    // fOut.AddOperation(op, {{"accuracy", std::to_string(0.001)}});

    // START OF NEW CODE
    adios2::Variable<double> tOut = bpIO.DefineVariable<double>("time");
    adios2::Variable<double> deltaTOut = bpIO.DefineVariable<double>("deltaT");
    adios2::Variable<double> hOut = bpIO.DefineVariable<double>("h");
    adios2::Variable<int> StepMaxOut = bpIO.DefineVariable<int>("MaxStep");
    adios2::Variable<int> vStep = bpIO.DefineVariable<int>("step");

    const std::string extent = "0 " + std::to_string(leny - 1) + " 0 " + std::to_string(lenz - 1) + " 0 " + std::to_string(lenx - 1);

    const std::string imageData = R"(

<?xml version="1.0"?>

<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">

<ImageData WholeExtent=")" + extent +
                                  R"(" Origin="0 0 0" Spacing="1 1 1">
<Piece Extent=")" + extent +
                                  R"(">
<PointData Scalars="F">

<DataArray Name="F" />

<DataArray Name="Laplace" />

<DataArray Name="TIME"> step 
</DataArray>

</PointData>

</Piece>

</ImageData>

</VTKFile>)";

    bpIO.DefineAttribute<std::string>("vtk.xml", imageData);
    bpIO.DefineAttribute<std::string>("meow", "meow meow ");

    // END OF NEW CODE
    gencore(rank, h, start_y, x, z, y); // gencore fills in the values of the x, z, y array on each rank

    for (int t = 1; t <= nsteps; t++)
    {
        double time = t * dt;
        evalF(num, rank, t, x, z, y, F); // evalF fills in the value of the two dimensional array F

        bpWriter.BeginStep();
        /** Put variables for buffering, template type is optional */
        double start = MPI_Wtime();
        bpWriter.Put(xOut, x.data());
        bpWriter.Put(zOut, z.data());
        bpWriter.Put(yOut, y.data());
        bpWriter.Put(tOut, time);
        bpWriter.Put(fOut, F.data());
        bpWriter.Put(vStep, t);
        if (t == 1)
        {
            bpWriter.Put(deltaTOut, dt);
            bpWriter.Put(hOut, h);
            bpWriter.Put(StepMaxOut, (int)nsteps);
        }
        bpWriter.EndStep();
        double stop = MPI_Wtime();
        double mpitime = stop - start;
        if (rank == 0)
            cout << "elapsed time = " << mpitime << endl;
    }

    /** Create bp file, engine becomes unreachable after this*/
    bpWriter.Close();

    if (rank == 0)

    {
        std::cout << "Ethan, I wrote file " << filename
                  << " to disk.\n";
    }

#if ADIOS2_USE_MPI

    MPI_Finalize();

#endif

    return 0;
}
