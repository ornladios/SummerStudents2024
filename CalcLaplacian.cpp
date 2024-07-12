#include <adios2.h>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <math.h>
#include <numeric>
#include <stdexcept> //std::invalid_argument std::exception
#include <stdio.h>
#include <vector>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif
#include <mpi.h>
using namespace std;

void calcLapace(int rank, int size, int nx, int nz, int ny, double h, std::vector<double> &F,
                std::vector<double> &laplace)
{
    // the whole ordering needs to be chnaged
    // it starts with nx,ny,nz
    // I THINK ALL THAT NEEDS TO BE CHNAGED IS NZ TO NY
    // BUT CHECK!!!!!!!!!!!!!!!!!
    // L = i+j*nx+k*nx*ny
    // l = nx*ny
    MPI_Request request;
    MPI_Status status;
    MPI_Request req_send_forward;
    MPI_Request req_send_backward;
    std::vector<double> forward_neighbor(nx * nz);  // sending forward
    std::vector<double> backward_neighbor(nx * nz); // sending backward
    forward_neighbor.resize(nx * nz, 0);
    backward_neighbor.resize(nx * nz, 0);
    std::vector<double> sending_bufferF; // buffer
    std::vector<double> sending_bufferB; // buffer
    // ny = rows
    // nx = cols
    // nz = width
    // THE ABOVE IS OLD LOOK AT THE FIRST COMMENT IN THE METHOD
    for (int k = 0; k < nz; k++) // loop needs to be changed
    {
        for (int i = 0; i < nx; i++)
        {
            int j = ny - 1;                   // change to z
            int L = i + nx * k + j * nx * nz; // global movemnent in cube
            int l = i + nx * k;               // local movement in plane
            forward_neighbor[l] = F[L];
        }
    }
    for (int k = 0; k < nz; k++) // loop needs to be changed
    {
        for (int i = 0; i < nx; i++)
        {
            int j = 0;                        // chnage to z
            int L = i + nx * k + j * nx * nz; // global movemnent in cube
            int l = i + nx * k;               // local movement in plane
            backward_neighbor[l] = F[L];
        }
    }
    sending_bufferF = forward_neighbor;
    sending_bufferB = backward_neighbor;
    if (size > 0)
    {
        if (rank == 0)
        { // chnage nx*nz to nx*ny for all

            MPI_Isend(sending_bufferF.data(), nx * nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD,
                      &req_send_forward);
            MPI_Recv(forward_neighbor.data(), nx * nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD,
                     &status);
        }
        else if (rank == size - 1)
        {
            sending_bufferB = backward_neighbor;
            MPI_Isend(sending_bufferB.data(), nx * nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                      &req_send_backward);
            MPI_Recv(backward_neighbor.data(), nx * nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                     &status);
        }
        else
        {
            MPI_Isend(sending_bufferB.data(), nz * nx, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                      &req_send_backward);
            MPI_Isend(sending_bufferF.data(), nx * nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD,
                      &req_send_forward);
            MPI_Recv(backward_neighbor.data(), nx * nz, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(forward_neighbor.data(), nx * nz, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD,
                     &status);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // calculates the laplace equation in the interior FOR ALL OF THEM!!!!!!!!
    // change order
    for (int j = 1; j < ny - 1; j++) // ny = row index =j
    {
        for (int k = 1; k < nz - 1; k++) // nz = width =k
        {
            for (int i = 1; i < nx - 1; i++) // nx = col =i
            {
                int L = i + nx * k + j * nx * nz;
                int ljp1 = i + nx * k + (j + 1) * nx * nz;
                int ljm1 = i + nx * k + (j - 1) * nx * nz;
                int lip1 = (i + 1) + nx * k + j * nx * nz;
                int lim1 = (i - 1) + nx * k + j * nx * nz;
                int lkp1 = i + nx * (k + 1) + j * nx * nz;
                int lkm1 = i + nx * (k - 1) + j * nx * nz;
                laplace[L] =
                    (F[lip1] + F[lim1] + F[ljp1] + F[ljm1] + F[lkp1] + F[lkm1] - 6 * F[L]) /
                    (h * h);
            }
        }
    }
    if (rank == 0)
    {
        // change order
        for (int k = 1; k < nz - 1; k++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                int j = ny - 1; // row index
                int l = i + nx * k;
                int L = i + nx * k + j * nx * nz;          // cols *row_index + col = location
                int ljp1 = i + nx * k + (j + 1) * nx * nz; // row + 1
                int ljm1 = i + nx * k + (j - 1) * nx * nz; // row - 1
                int lip1 = (i + 1) + nx * k + j * nx * nz; // column + 1
                int lim1 = (i - 1) + nx * k + j * nx * nz; // column - 1
                int lkp1 = i + nx * (k + 1) + j * nx * nz;
                int lkm1 = i + nx * (k - 1) + j * nx * nz;

                laplace[L] = (F[lip1] + F[lim1] + forward_neighbor[l] + F[ljm1] + F[lkp1] +
                              F[lkm1] - 6 * F[L]) /
                             (h * h);
            }
        }
    }
    else if (rank == size - 1) // last proc
    {
        // chnage order
        for (int k = 1; k < nz - 1; k++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                int j = 0; // row index
                int l = i + nx * k;
                int L = i + nx * k + j * nx * nz;          // cols *row_index + col = location
                int ljp1 = i + nx * k + (j + 1) * nx * nz; // row + 1
                int ljm1 = i + nx * k + (j - 1) * nx * nz; // row - 1
                int lip1 = (i + 1) + nx * k + j * nx * nz; // column + 1
                int lim1 = (i - 1) + nx * k + j * nx * nz; // column - 1
                int lkp1 = i + nx * (k + 1) + j * nx * nz;
                int lkm1 = i + nx * (k - 1) + j * nx * nz;

                laplace[L] = (F[lip1] + F[lim1] + F[ljp1] + backward_neighbor[l] + F[lkp1] +
                              F[lkm1] - 6 * F[L]) /
                             (h * h);
            }
        }
    }
    else
    { // chnage order
        for (int k = 1; k < nz - 1; k++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                int j = 0; // row index
                int l = i + nx * k;
                int L = i + nx * k + j * nx * nz;          // cols *row_index + col = location
                int ljp1 = i + nx * k + (j + 1) * nx * nz; // row + 1
                int ljm1 = i + nx * k + (j - 1) * nx * nz; // row - 1
                int lip1 = (i + 1) + nx * k + j * nx * nz; // column + 1
                int lim1 = (i - 1) + nx * k + j * nx * nz; // column - 1
                int lkp1 = i + nx * (k + 1) + j * nx * nz;
                int lkm1 = i + nx * (k - 1) + j * nx * nz;

                laplace[L] = (F[lip1] + F[lim1] + F[ljp1] + backward_neighbor[l] + F[lkp1] +
                              F[lkm1] - 6 * F[L]) /
                             (h * h);
            }
        }
        // chang order
        for (int k = 1; k < nz - 1; k++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                int j = ny - 1; // row index
                int l = i + nx * k;
                int L = i + nx * k + j * nx * nz;          // cols *row_index + col = location
                int ljp1 = i + nx * k + (j + 1) * nx * nz; // row + 1
                int ljm1 = i + nx * k + (j - 1) * nx * nz; // row - 1
                int lip1 = (i + 1) + nx * k + j * nx * nz; // column + 1
                int lim1 = (i - 1) + nx * k + j * nx * nz; // column - 1
                int lkp1 = i + nx * (k + 1) + j * nx * nz;
                int lkm1 = i + nx * (k - 1) + j * nx * nz;

                laplace[L] = (F[lip1] + F[lim1] + forward_neighbor[l] + F[ljm1] + F[lkp1] +
                              F[lkm1] - 6 * F[L]) /
                             (h * h);
            }
        }
    }
}

int main(int argc, char **argv)
{
    int rank, size;
#if ADIOS2_USE_MPI
    // std::cout << "Enter file name, y size, z size, x size, time steps, function number 1-4" <<
    // endl;
    int provided;
    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    // rank is giving a number to each proccess
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //  the size of something
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#else
// std::cout << "Enter file name, y size, z size, x size, time steps, function number 1-4" <<
// endl;
std:
    cerr << "Error: This code is designed to compiled with mpi!" << std::endl;
    return -1;
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD);
#else
    adios2::ADIOS adios("adios2.xml");
#endif

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <inputFname> <outputFilename>" << endl;
        return -1;
    }
    std::string inputFname = argv[1];
    std::string outputFilename = argv[2];

    auto inIO = adios.DeclareIO("WriteIO");
    auto reader = inIO.Open(inputFname, adios2::Mode::Read);

    // Creating a new IO and engine for the output file
    adios2::IO bpIO = adios.DeclareIO("LaplaceOutput");
    adios2::Engine bpWriter = bpIO.Open(outputFilename, adios2::Mode::Write);

    // Defining the variable lOut, Maxstep, delta T, x,y,z to write in
    adios2::Variable<double> lOut;
    adios2::Variable<int> StepMaxOut;
    adios2::Variable<double> deltaTOut;
    adios2::Variable<double> xOut;
    adios2::Variable<double> yOut;
    adios2::Variable<double> zOut;

    // Defining the variable h
    double h;
    int step = 0;
    // these variables are needed for the time derivative
    double deltaT;
    int MaxStep;

    while (true)
    {

        // Get the status of the reader
        auto status = reader.BeginStep();

        // If status is not ok (such as when there is nothing left to be read), break out of the while loop
        if (status != adios2::StepStatus::OK)
        {
            break;
        }

        // Inquire the variable F from the input IO
        auto varF = inIO.InquireVariable<double>("F");
        auto varX = inIO.InquireVariable<double>("x");
        auto varY = inIO.InquireVariable<double>("y");
        auto varZ = inIO.InquireVariable<double>("z");
        // If it is the first step, get the value of h, since it stays constant
        if (step == 0)
        {
            auto varh = inIO.InquireVariable<double>("h");
            reader.Get<double>(varh, &h, adios2::Mode::Sync);
            auto varT = inIO.InquireVariable<double>("deltaT");
            reader.Get<double>(varT, &deltaT, adios2::Mode::Sync);
            auto varMaxStep = inIO.InquireVariable<int>("MaxStep");
            reader.Get<int>(varMaxStep, &MaxStep, adios2::Mode::Sync);

            std::cout << "Max step: " << MaxStep << endl;
            std::cout << "Delta T: " << deltaT << endl;
        }

        // Getting the shape of F (which should in the format of y, x, z) NOTE: change accordingly
        // to coordinate system

        auto shapeF = varF.Shape();
        size_t lenx = shapeF[0];
        size_t leny = shapeF[1]; // y is the size that is changing
        size_t lenz = shapeF[2];

        auto shapeX = varX.Shape();
        auto shapeY = varY.Shape();
        auto shapeZ = varZ.Shape();

        std::vector<double> arrX(lenx, 0.0);
        std::vector<double> arrY(leny, 0.0);
        std::vector<double> arrZ(lenz, 0.0);

        varX.SetSelection({{0}, {lenx}});
        reader.Get<double>(varX, arrX.data(), adios2::Mode::Sync);

        varY.SetSelection({{0}, {leny}});
        reader.Get<double>(varY, arrY.data(), adios2::Mode::Sync);

        varZ.SetSelection({{0}, {lenz}});
        reader.Get<double>(varZ, arrZ.data(), adios2::Mode::Sync);

        // Calculating the local_len_y and the remainder
        // NOTE: This is assuming that the y dimension is the dimension that is being split!!!
        size_t local_len_y = leny / size;
        size_t remainder = leny % size;
        if (rank < remainder)
        {
            local_len_y++;
        }

        // Calculating the total length of the global and local arrays
        size_t len_global = lenx * leny * lenz;
        size_t len_local = local_len_y * lenx * lenz;

        // Calculating the start_y index of the local array
        size_t start_y = rank * (leny / size) + (rank < remainder ? rank : remainder);

        // Defining the variable lOut if it is the first step
        if (step == 0)
        {
            lOut = bpIO.DefineVariable<double>("Laplace", {lenx, leny, lenz}, {0, start_y, 0},
                                               {lenx, local_len_y, lenz}, adios2::ConstantDims);
            StepMaxOut = bpIO.DefineVariable<int>("MaxStep");
            deltaTOut = bpIO.DefineVariable<double>("deltaT");
            xOut = bpIO.DefineVariable<double>("x", {lenx}, {0}, {lenx}, adios2::ConstantDims);
            yOut = bpIO.DefineVariable<double>("y", {leny}, {0}, {leny}, adios2::ConstantDims); // weird former changing size
            zOut = bpIO.DefineVariable<double>("z", {lenz}, {0}, {lenz}, adios2::ConstantDims);
        }

        // Creating vectors of size len_local to store the incoming data
        std::vector<double> arrF(len_local, 0.0);
        std::vector<double> laplace(len_local, 0.0);

        // Selecting which section of the data to read
        varF.SetSelection({{start_y, 0, 0}, {local_len_y, lenz, lenx}});

        // Reading the data from the bpFile into arrF
        reader.Get<double>(varF, arrF.data(), adios2::Mode::Sync);

        // Calculating the laplace equation
        // NOTE: CURRENTLY ONLY WORKS WITH MORE THAN 1 PROCESSORS (COULD BE CHANGED TO WORK WITH 1 PROCESSOR...TECHNICALLY)
        if (size > 1)
        {
            calcLapace(rank, size, lenx, lenz, local_len_y, h, arrF, laplace);
        }
        else
        {
            std::cout << "Error: This code is designed to run with multiple processors" << std::endl;
            break;
        }

        reader.EndStep();

        bpWriter.BeginStep();
        if (step == 0)
        {

            bpWriter.Put(deltaTOut, deltaT);   // something wrong
            bpWriter.Put(StepMaxOut, MaxStep); // something wrong
        }

        // Writing out the Laplace data to outputFilename

        double start = MPI_Wtime();
        bpWriter.Put(lOut, laplace.data());
        bpWriter.Put(xOut, arrX.data());
        bpWriter.Put(yOut, arrY.data());
        bpWriter.Put(zOut, arrZ.data());
        bpWriter.EndStep();
        double stop = MPI_Wtime();
        double mpitime = stop - start;

        if (rank == 0)
            cout << "elapsed time = " << mpitime << " step: " << step << endl;

        step++;
    }

    reader.Close();
    bpWriter.Close();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}