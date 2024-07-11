/*
 * Notes:
 *
 * In this code, slices in a 3D grid, is similar to a loaf of bread position vertically in front of
 * you, and you are cutting horizontally to create slices (should result in normal bread slices).
 *
 * The slice closest to you will have an index of 0.
 *
 * The width of the slice is the x-axis, the height of the slice is the y-axis, and the depth of the
 * slice is the z-axis.
 *
 *
 * Here's a 2d visualization (I tried my best):

 *  The coordinate system is a bit weird in this format its (z, x, y) instead of (x, y, z).
*/

#include <adios2.h>
#include <ios>       //std::ios_base::failure
#include <iostream>  //std::cout
#include <math.h>    // This allows me to use cosine
#include <stdexcept> //std::invalid_argument std::exception
#include <stdio.h>   // This allows me to use printf
#include <vector>
#if ADIOS2_USE_MPI

#include <mpi.h>

#endif
#define PI 3.141592653589793
using namespace std;

void genCoords(int rank, int size, double start_x, double end_x, std::vector<double> &xvector,
               std::vector<double> &yvector, std::vector<double> &zvector, int start_z, int len_z)
{
    if (start_x >= end_x)
    {
        throw std::invalid_argument("start_x must be less than end_x");
    }

    // z-axis (depth) grows as processor increases -> we need to scale zvector by rank
    double delta_z = (end_x - start_x) / (len_z - 1);
    for (int i = 0; i < zvector.size(); i++)
        zvector[i] = (start_z + i) * delta_z;

    // doesnt grow as process size increases
    double delta_y = (end_x - start_x) / (yvector.size() - 1);
    for (int i = 0; i < yvector.size(); i++)
        yvector[i] = i * delta_y;

    // doesnt grow as process size increases
    double delta_x = (end_x - start_x) / (xvector.size() - 1);
    for (int i = 0; i < xvector.size(); i++)
        xvector[i] = i * delta_x;
}

void calcF(int num, int rank, double t, std::vector<double> &x, std::vector<double> &z,
           std::vector<double> &y, std::vector<double> &F)
{
    double theta = fmod(t, 2.0 * PI);
    double phi = fmod(t / 2.0, 2.0 * PI); // for nums
    double x0 = sin(phi) + 3;
    double y0 = sin(phi) + 3;
    double z0 = PI;
    double spin = -1;
    double x1 = spin * sin(4 * phi) + x0;
    double y1 = spin * cos(4 * phi) + y0;
    double z1 = PI;

    for (int k = 0; k < z.size(); k++)
    {
        for (int i = 0; i < x.size(); i++) // i =col
        {

            for (int j = 0; j < y.size(); j++) //  j= row
            {

                int L = i + x.size() * k + j * x.size() * z.size();
                if (num == 1) // boring
                {
                    F[L] = x[i] * x[i] + y[j] * y[j] + z[k] * z[k];
                }
                else if (num == 2) // cool
                {
                    F[L] = 1.0 / (sqrt(pow(x[i] - x0, 2) + pow(y[j] - y0, 2) + pow(z[k] - z0, 2)));
                }
                else if (num == 3) // intresting
                {
                    F[L] = 1.0 / (sqrt(pow(x[i] - x0, 2) + pow(y[j] - y0, 2) + pow(z[k] - z0, 2))) +
                           1.0 / (sqrt(pow(x[i] - x1, 2) + pow(y[j] - y1, 2) + pow(z[k] - z1, 2)));
                }
                else // WOWERS
                {
                    F[L] = 1.5 / (pow((x[i] - x0) * cos(theta) - (y[j] - y0) * sin(theta), 2) +
                                  4 * pow((x[i] - x0) * sin(theta) + (y[j] - y0) * cos(theta), 2) +
                                  pow(z[k] - z0, 2)) +
                           0.5 / (pow((x[i] - x1) * cos(theta) - (y[j] - y1) * sin(theta), 2) +
                                  4 * pow((x[i] - x1) * sin(theta) + (y[j] - y1) * cos(theta), 2) +
                                  pow(z[k] - z1, 2));
                }
            }
        }
    }
}
int main(int argc, char *argv[])
{
    int rank, size;
#if ADIOS2_USE_MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#else
    rank = 0;
    size = 1;

#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios("./adios2.xml", MPI_COMM_WORLD);
#else
    adios2::ADIOS adios("./adios2.xml");
#endif
    // Making sure our user passed in the correct number of arguments
    // We need 6 arguments: filename, length_z, length_x, length_y, tmax
    // We check for 7 because argv[0] is always the name of the program
    if (argc < 7)
    {
        if (rank == 0)
        {
            // Only rank 0 prints the error message to avoid duplication
            fprintf(stderr, "Usage: %s <filename> <length_z> <length_x> <length_y> <tmax> <num>\n",
                    argv[0]);
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1); // Abort the program with a non-zero exit code
#endif
        return -1;
    }

    // Parsing the cmdline arguments into variables
    std::string filename;
    filename = argv[1];           // name of the .bp file output
    size_t len_z = atoi(argv[2]); // global length
    size_t len_x = atoi(argv[3]); // global length
    size_t len_y = atoi(argv[4]); // global length
    size_t tmax = atoi(argv[5]);
    if (!(len_z == len_x && len_z == len_y && len_x == len_y))
    {
        if (rank == 0)
        {
            // Only rank 0 prints the error message to avoid duplication
            fprintf(stderr, "Error: len_z, len_x, and len_y must be equal to each other.\n");
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1); // Abort the program with a non-zero exit code
#endif
        return -1;
    }

    int num = atoi(argv[6]); // num is choose function

    // Total number of values in the 2d grid
    size_t local_len_z = len_z / size;
    size_t remainder = len_z % size;

    // If the rank is less than the remainder, we need to add one more to the local_len_z, to
    // handle global length that can't be evenly divided by the number of processes
    //
    // Ex. length of 100 among 3 processors
    // First one will have 34
    // Second one will have 33
    // Third one will have 33

    if (rank < remainder)
    {
        local_len_z++;
    }

    size_t len_global = len_x * len_y * len_z;
    size_t len_local = local_len_z * len_x * len_y;

    // starting z index for this rank
    size_t start_z = rank * (len_z / size) + (rank < remainder ? rank : remainder);

    // Instantiate the x and y vectors
    std::vector<double> x, y, z;

    // Initilize the x, y, and z vectors with their respective lengths
    x.resize(len_x, 0);
    y.resize(len_y, 0);
    z.resize(local_len_z, 0);

    // Starting and end-point of the x values
    double start_x = 0.0;
    double end_x = 2 * PI;

    // Initialize a double vector of size len_global
    std::vector<double> F(len_local, 0.0); // F is the function we are trying to differentiate

    // Creating the adios2 IO object
    adios2::IO bpIO = adios.DeclareIO("WriteIO");

    // Opening <filename> for writing with adios2 bpIO object
    adios2::Engine bpWriter = bpIO.Open(filename, adios2::Mode::Write);
    adios2::Operator op = adios.DefineOperator("mgard", "mgard");

    // Defining a variable for the x values
    adios2::Variable<double> xOut =
        bpIO.DefineVariable<double>("x", {len_x}, {0}, {len_x}, adios2::ConstantDims);

    // Defining a variable for the y values
    adios2::Variable<double> yOut =
        bpIO.DefineVariable<double>("y", {len_y}, {0}, {len_y}, adios2::ConstantDims);

    // Defining a variable for the z values
    adios2::Variable<double> zOut =
        bpIO.DefineVariable<double>("z", {len_z}, {start_z}, {local_len_z}, adios2::ConstantDims);

    adios2::Variable<double> tOut = bpIO.DefineVariable<double>("time");
    adios2::Variable<double> deltaTOut = bpIO.DefineVariable<double>("deltaT");
    adios2::Variable<double> hOut = bpIO.DefineVariable<double>("h");
    adios2::Variable<int> tmaxOut = bpIO.DefineVariable<int>("MaxStep");
    adios2::Variable<int> stepOut = bpIO.DefineVariable<int>("step");

    // Defining a variable for the function values
    // "F" is going to have a global size of (len_z * size) by len_x by len_z
    // Each process is gonna start writing at index [rank * len_z][0][0] meaning rank * len_z
    // slice of the 3d grid Each process is gonna write (len_z, len_x, len_y) values
    adios2::Variable<double> fOut =
        bpIO.DefineVariable<double>("F", {len_z, len_x, len_y}, {start_z, 0, 0},
                                    {local_len_z, len_x, len_y}, adios2::ConstantDims);

    const std::string extent = "0 " + std::to_string(size * z.size() - 1) + " 0 " +
                               std::to_string(x.size() - 1) + " 0 " + std::to_string(y.size() - 1);

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

    // Change in value between points in all axis
    double dt = 1.0 / tmax;
    double h = (end_x - start_x) / (len_x - 1);

    // Generate the x and y coordinates
    genCoords(rank, size, start_x, end_x, x, y, z, start_z, len_z);

    for (int t = 1; t <= tmax; t++)
    {
        double time = t * dt;
        // Calculate the function values
        calcF(num, rank, t, x, z, y, F);

#if ADIOS2_USE_MPI
        double start = MPI_Wtime();
#endif
        // Begin the step
        bpWriter.BeginStep();

        // Put variables
        if (rank == 0)
        {
            bpWriter.Put(xOut, x.data()); // writing the x values
            bpWriter.Put(yOut, y.data()); // writing the y values
            bpWriter.Put(zOut, z.data()); // writing the z values
            bpWriter.Put(tOut, time);     // writing out the time used to calculate F
            bpWriter.Put(stepOut, t);
        }

        bpWriter.Put(fOut, F.data()); // this must be written for each rank

        if (t == 1)
        {
            bpWriter.Put(deltaTOut, dt);
            bpWriter.Put(hOut, h);
            bpWriter.Put(tmaxOut, (int)tmax);
        }

        // End the step
        bpWriter.EndStep();

#if ADIOS2_USE_MPI
        double end = MPI_Wtime();
        if (rank == 0)
        {
            std::cout << "Time Elapsed: " << end - start << " seconds" << std::endl;
        }
#endif
    }

    // Close the writer
    bpWriter.Close();

    // Print a message to the user, that all processes have finished writing
    if (rank == 0)
    {
        std::cout << "Wrote file " << filename
                  << " to disk. It can now be read by running "
                     "`bpls -l "
                  << filename << "`.\n";
    }

#if ADIOS2_USE_MPI

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

#endif

    return 0;
}