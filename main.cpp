#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>

using OptionMap = boost::program_options::variables_map;

auto getOptions(int argc, char *argv[]) -> OptionMap
{
  namespace po = boost::program_options;

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")("mesh,m", po::value<std::string>(), "The partitioned mesh prefix used as input (only VTU format is accepted)(Looking for <prefix>_<#filerank>.vtu).")("output,o", po::value<std::string>(), "The output mesh. Can be VTK or VTU format. If it is not given <inputmesh>_joined.vtk will be used.")("recovery,r", po::value<std::string>(), "The path to the recovery file to fully recover it's state.")("numparts,n", po::value<size_t>()->default_value(0), "The number of parts to read from the input mesh. By default the entire mesh is read.")("directory,dir", po::value<std::string>()->default_value("."), "Directory for output files (optional)");

  po::variables_map vm;
  try {
    po::store(parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      std::exit(EXIT_SUCCESS);
    }
    // Needs to be called
    po::notify(vm);

    if (!vm.count("mesh")) {
      std::cout << "You must specify a mesh file to read from." << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    std::cerr << desc << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return vm;
}

auto count_partitions(const std::string &prefix) -> size_t
{
  namespace fs = boost::filesystem;

  std::string filename;
  size_t      count = 0;

  while (true) {
    filename = prefix + "_" + std::to_string(count) + ".vtu";
    if (!fs::exists(filename)) {
      break;
    }
    ++count;
  }
  return count;
}

void write_mesh(const std::string &filename, vtkSmartPointer<vtkUnstructuredGrid> mesh)
{
  auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(mesh);
  writer->Write();
}

auto partitionwise_merge(const std::string &prefix, size_t numparts) -> vtkSmartPointer<vtkUnstructuredGrid>
{
  auto             joined_mesh   = vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto             joined_points = vtkSmartPointer<vtkPoints>::New();
  auto             joined_cells  = vtkSmartPointer<vtkCellArray>::New();
  std::vector<int> joined_cell_types;

  std::vector<vtkSmartPointer<vtkDoubleArray>> joined_data_vec;
  std::vector<std::string>                     joined_data_names;

  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  for (size_t i = 0; i < numparts; ++i) {
    // Read mesh
    auto partname = prefix + "_" + std::to_string(i) + ".vtu";
    reader->SetFileName(partname.c_str());
    reader->Update();
    // Extract mesh
    auto grid = reader->GetOutput();
    // Cells
    const auto offset    = joined_points->GetNumberOfPoints();
    auto       num_cells = grid->GetCells()->GetNumberOfCells();
    joined_cell_types.reserve(joined_cell_types.size() + num_cells);
    auto cellIds = vtkSmartPointer<vtkIdList>::New();
    for (int j = 0; j < num_cells; ++j) {
      cellIds->Reset();
      std::for_each(grid->GetCell(j)->GetPointIds()->begin(), grid->GetCell(j)->GetPointIds()->end(), [&cellIds, &offset](auto &pointId) { cellIds->InsertNextId(pointId + offset); });
      joined_cells->InsertNextCell(cellIds);
      joined_cell_types.push_back(grid->GetCellType(j));
    }
    // Points
    auto points = grid->GetPoints();
    joined_points->InsertPoints(joined_points->GetNumberOfPoints(), grid->GetNumberOfPoints(), 0, points);
    //  Point Data
    auto part_point_data = grid->GetPointData();
    auto num_arrays      = part_point_data->GetNumberOfArrays();
    for (int j = 0; j < num_arrays; ++j) {
      auto part_data = part_point_data->GetArray(j);
      auto name      = part_data->GetName();
      if (std::find(joined_data_names.begin(), joined_data_names.end(), name) == joined_data_names.end()) {
        joined_data_names.emplace_back(name);
        auto new_joined_data = vtkSmartPointer<vtkDoubleArray>::New();
        new_joined_data->SetName(name);
        new_joined_data->SetNumberOfComponents(part_data->GetNumberOfComponents());
        joined_data_vec.push_back(new_joined_data);
      }
      auto joined_data = joined_data_vec[j];
      joined_data->InsertTuples(joined_data->GetNumberOfTuples(), part_data->GetNumberOfTuples(), 0, part_data);
    }
  }

  joined_mesh->SetPoints(joined_points);
  for (const auto &data : joined_data_vec) {
    joined_mesh->GetPointData()->AddArray(data);
  }
  joined_mesh->SetCells(joined_cell_types.data(), joined_cells);

  return joined_mesh;
}

void join(int argc, char *argv[])
{
  namespace fs        = boost::filesystem;
  auto        options = getOptions(argc, argv);
  std::string prefix  = options["mesh"].as<std::string>();
  std::string output{};
  std::string recovery{};
  if (options.find("output") == options.end()) {
    output = prefix + "_joined.vtk";
  } else {
    output = options["output"].as<std::string>();
  }
  if (options.find("recovery") != options.end()) {
    recovery = options["recovery"].as<std::string>();
  } else {
    recovery = prefix + "_recovery.json";
  }
  std::string directory = options["directory"].as<std::string>();
  size_t      numparts  = options["numparts"].as<size_t>();

  if (numparts == 0) {
    numparts = count_partitions(prefix);
  }

  if (fs::exists(recovery)) {
    std::cout << "Recovery file found. Will try to recover the state." << std::endl;
  } else {
    std::cout << "Recovery file not found. Partition-wise merging will be done." << std::endl;
  }

  auto joined_mesh = partitionwise_merge(prefix, numparts);
  write_mesh(output, joined_mesh);
}

void recovery_merge()
{
}

auto main(int argc, char *argv[]) -> int
{
  join(argc, argv);
  return EXIT_SUCCESS;
}