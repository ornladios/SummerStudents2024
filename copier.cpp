#include <adios2.h>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <map>
#include <math.h>
#include <mpi.h>
#include <numeric>
#include <stdexcept> //std::invalid_argument std::exception
#include <stdio.h>
#include <vector>

int main(int argc, char **argv)
{
    int rank = 0, numProcs = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (argc != 3)
    {
        std::cerr << "****************************************************" << std::endl;
        std::cerr << "Error: Usage: " << argv[0] << " input output" << std::endl;
        std::cerr << "  Please rerun with the proper arguments." << std::endl;
        std::cerr << "  Thanks for stopping by!" << std::endl << std::endl;
        MPI_Finalize();
        return 0;
    }
    if (numProcs != 1)
    {
        std::cerr << "****************************************************" << std::endl;
        std::cerr << " Error. This probably only works with 1 rank..." << std::endl << std::endl;
        MPI_Finalize();
        return 0;
    }

    adios2::ADIOS adios("adios2.xml", MPI_COMM_WORLD);

    std::string inputFname(argv[1]), outputFname(argv[2]);

    if (rank == 0)
        std::cerr << "Reading from: " << inputFname << std::endl;

    auto inIO = adios.DeclareIO("WriteIO");
    auto reader = inIO.Open(inputFname, adios2::Mode::Read);

    auto outIO = adios.DeclareIO("Output");
    auto writer = outIO.Open(outputFname, adios2::Mode::Write);

    int step = 0;
    while (true)
    {
        std::vector<
            std::pair<std::string, std::pair<std::vector<double>, std::vector<std::size_t>>>>
            data;
        std::vector<std::pair<std::string, std::pair<std::vector<int>, std::vector<std::size_t>>>>
            datai;
        std::map<std::string, std::string> attrs;

        auto status = reader.BeginStep();
        if (status != adios2::StepStatus::OK)
            break;

        auto variables = inIO.AvailableVariables();
        std::cout << "Reading Step= " << step << std::endl;
        for (const auto &vi : variables)
        {
            if (vi.first == "step")
            {
                auto var = inIO.InquireVariable<int>(vi.first);
                if (!var)
                    continue;
                auto shape = var.Shape();
                std::size_t dataSz =
                    std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<std::size_t>());
                std::vector<int> arr(dataSz, 0.0);
                reader.Get<int>(var, arr.data(), adios2::Mode::Sync);
                std::vector<std::size_t> varS = shape;
                datai.push_back(std::move(std::make_pair(var.Name(), std::make_pair(arr, varS))));
            }
            else
            {
                auto var = inIO.InquireVariable<double>(vi.first);
                if (!var)
                    continue;

                auto shape = var.Shape();
                std::size_t dataSz =
                    std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<std::size_t>());
                // std::cout<<"  Read "<<vi.first<<" : sz= "<<dataSz<<std::endl;
                std::vector<double> arr(dataSz, 0.0);
                reader.Get<double>(var, arr.data(), adios2::Mode::Sync);
                std::vector<std::size_t> varS = shape;
                data.push_back(std::move(std::make_pair(var.Name(), std::make_pair(arr, varS))));
            }
        }

        if (step == 0)
        {
            auto allAttrs = inIO.AvailableAttributes();
            for (const auto &ai : allAttrs)
            {
                const auto a = inIO.InquireAttribute<std::string>(ai.first);
                attrs[a.Name()] = a.Data()[0];
                // std::cout<<" attr: "<<a.Name()<<" "<<a.Data()[0]<<std::endl;
            }
        }
        reader.EndStep();

        // Write data.
        std::cout << "Writing Step= " << step << std::endl;
        status = writer.BeginStep();
        if (status != adios2::StepStatus::OK)
        {
            std::cerr << "Failure on writer.BeginStep()" << std::endl;
            break;
        }

        if (step == 0 && !attrs.empty())
        {
            for (const auto &a : attrs)
            {
                // std::cout<<"      attr= "<<a.first<<std::endl;
                outIO.DefineAttribute<std::string>(a.first, a.second);
            }
        }

        for (const auto &d : data)
        {
            auto varNm = d.first;
            auto varSz = d.second.second;
            std::size_t sz = d.second.first.size();
            // const adios2::Dims shape, start, count; //{sz}, start{0}, count{sz};
            std::vector<std::size_t> shape, start, count;
            for (std::size_t i = 0; i < varSz.size(); i++)
            {
                shape.push_back(varSz[i]);
                count.push_back(varSz[i]);
                start.push_back(0);
            }
            // const adios2::Box<adios2::Dims> sel(0,sz); //(start, count);
            // const adios2::Box<adios2::Dims> sel(start, count); //(start, count);
            adios2::Variable<double> vOut = outIO.InquireVariable<double>(varNm);
            if (!vOut)
                vOut = outIO.DefineVariable<double>(varNm, shape, start, count);

            // vOut.SetSelection(sel);
            // std::cout<<"  Write "<<varNm<<" : sz= "<<sz<<std::endl;
            writer.Put(vOut, d.second.first.data());
        }
        for (const auto &d : datai)
        {
            auto varNm = d.first;
            auto varSz = d.second.second;
            std::size_t sz = d.second.first.size();
            // const adios2::Dims shape, start, count; //{sz}, start{0}, count{sz};
            std::vector<std::size_t> shape, start, count;
            for (std::size_t i = 0; i < varSz.size(); i++)
            {
                shape.push_back(varSz[i]);
                count.push_back(varSz[i]);
                start.push_back(0);
            }

            // const adios2::Box<adios2::Dims> sel(start, count); //(start, count);
            adios2::Variable<int> vOut = outIO.InquireVariable<int>(varNm);
            if (!vOut)
                vOut = outIO.DefineVariable<int>(varNm, shape, start, count);

            // vOut.SetSelection(sel);
            // std::cout<<"  Write "<<varNm<<" : sz= "<<sz<<std::endl;
            writer.Put(vOut, d.second.first.data());
        }
        writer.EndStep();
        step++;
    }

    reader.Close();
    writer.Close();

    MPI_Finalize();

    return 0;
}
