#include "Parameters.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

class OptionInterpreter {
 public:
  OptionInterpreter(const std::string &option_str, const std::string &value_str)
      : option_str_(option_str),
        value_str_(value_str) {
  }

  // Load a single value matching the option name
  template<typename T>
  bool load_value(const std::string &target_option_name,
                  T &target_option_value) const {
    if (option_str_ == target_option_name) {
      if (!load_value_impl(value_str_, target_option_value)) {
        std::cerr << "Error loading option: " << target_option_name
                  << std::endl;
        return false;
      }

      return true;
    } else {
      return false;
    }
  }

  // Load a sequence of values matching the option name
  template<typename T>
  bool load_values(const std::string &target_option_name,
                   std::vector<T> &target_option_values) const {
    if (option_str_ == target_option_name) {
      target_option_values.clear();
      std::istringstream istr(value_str_);
      T value;
      while (istr >> value) {
        target_option_values.push_back(value);
      }

      if (target_option_values.empty()) {
        std::cerr << "Error: no value loaded for option: " << target_option_name
                  << std::endl;
        return false;
      }
      return true;
    } else {
      return false;
    }
  }

  // Load an enum matching the option name
  template<typename EnumT>
  bool load_enum(const std::string &target_option_name, int enum_value_count,
                 EnumT &value) const {
    if (option_str_ == target_option_name) {
      int enum_int = 0;
      if (load_value_impl(value_str_, enum_int)) {
        if (enum_int >= 0 && enum_int < enum_value_count) {
          value = static_cast<EnumT>(enum_int);
          return true;
        }
      }

      std::cerr << "Error loading option: " << target_option_name << std::endl;
      return false;
    } else {
      return false;
    }
  }

 private:
  std::string option_str_, value_str_;

  bool load_value_impl(const std::string &str, double &value) const {
    try {
      value = std::stod(str);
    } catch (const std::invalid_argument& ia) {
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      return false;
    } catch (const std::out_of_range &oor) {
      std::cerr << "Out of Range error: " << oor.what() << std::endl;
      return false;
    }

    return true;
  }

  bool load_value_impl(const std::string &str, int &value) const {
    try {
      value = std::stoi(str);
    } catch (const std::invalid_argument& ia) {
      std::cerr << "Invalid argument: " << ia.what() << std::endl;
      return false;
    } catch (const std::out_of_range &oor) {
      std::cerr << "Out of Range error: " << oor.what() << std::endl;
      return false;
    }

    return true;
  }

  bool load_value_impl(const std::string &str, bool &value) const {
    int bool_value = 0;
    if (load_value_impl(str, bool_value)) {
      value = (bool_value != 0);
      return true;
    } else {
      return false;
    }
  }
};

// Load options from file
bool Parameters::load(const char* filename) {
  std::ifstream ifile(filename);
  if (!ifile.is_open()) {
    std::cerr << "Error while opening file " << filename << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(ifile, line)) {
    std::string::size_type pos = line.find_first_not_of(' ');
    if (pos == std::string::npos) {
      continue;
    }

    // Check for comment line
    else if (line.at(pos) == '#') {
      continue;
    }

    std::string::size_type first_pos = line.find_first_not_of(' ');
    std::string trimmed_line = line.substr(first_pos, std::string::npos);
    std::string::size_type end_pos = trimmed_line.find_first_of(' ');
    std::string option_str = trimmed_line.substr(pos, end_pos - pos);
    std::string value_str = trimmed_line.substr(end_pos + 1, std::string::npos);
    OptionInterpreter opt(option_str, value_str);

    if (!(opt.load_value("HeatSolverMaxIter", heat_solver_max_iter)
        || opt.load_value("HeatSolverEps", heat_solver_eps)
        || opt.load_value("HeatSolverConvergeCheckFrequency", heat_solver_convergence_check_frequency)
        || opt.load_value("GradSolverMaxIter", grad_solver_max_iter)
        || opt.load_value("GradSolverEps", grad_solver_eps)
        || opt.load_value("Penalty", penalty)
        || opt.load_value("GradSolverOutputFrequency", grad_solver_output_frequency)
        || opt.load_value("GradSolverConvergeCheckFrequency", grad_solver_convergence_check_frequency)
        || opt.load_values("SourceVertices", source_vertices))) {
      std::cerr << "Unable to parse option " << option_str << std::endl;
      return false;
    }
  }

  if(!valid_parameters()){
    std::cerr << "Error: invalid parameter value" << std::endl;
    return false;
  }

  std::sort(source_vertices.begin(), source_vertices.end());

  std::cout << "Successfully loaded options from file " << filename
            << std::endl;
  return true;
}

template<typename T>
bool check_lower_bound(const std::string &name, T val, T lower_bound,
                       bool allow_equal) {
  bool valid = allow_equal ? (val >= lower_bound) : (val > lower_bound);
  if (!valid) {
    std::cerr << "Error: " << name << " must be "
        << (allow_equal ? ("at least") : ("larger than")) << " " << lower_bound
        << std::endl;
  }

  return valid;
}

template<typename T>
bool check_upper_bound(const std::string &name, T val, T upper_bound,
                       bool allow_equal) {
  bool valid = allow_equal ? (val <= upper_bound) : (val < upper_bound);
  if (!valid) {
    std::cerr << "Error: " << name << " must be "
        << (allow_equal ? ("at most") : ("smaller than")) << " " << upper_bound
        << std::endl;
  }

  return valid;
}


bool check_nonempty_index_sequence(const std::string &name, const std::vector<int> &seq)
{
  bool valid = (!seq.empty());
  if(!valid){
    std::cerr << "Error: " << name << " is empty" << std::endl;
  }

  for(int i = 0; i < static_cast<int>(seq.size()); ++ i){
    if(seq[i] < 0){
      std::cerr << "Error: invalid index " << seq[i] << std::endl;
      valid = false;
      break;
    }
  }

  return valid;
}

// Check whether the parameter values are valid
bool Parameters::valid_parameters() const {
  return check_lower_bound("HeatSolverMaxIter", heat_solver_max_iter, 0, false)
      && check_lower_bound("HeatSolverEps", heat_solver_eps, 0.0, false)
      && check_lower_bound("HeatSolverConvergeCheckFrequency", heat_solver_convergence_check_frequency, 0, false)
      && check_lower_bound("GradSolverMaxIter", grad_solver_max_iter, 0, false)
      && check_lower_bound("GradSolverEps", grad_solver_eps, 0.0, false)
      && check_lower_bound("Penalty", penalty, 0.0, false)
      && check_lower_bound("GradSolverOutputFrequency", grad_solver_output_frequency, 0, false)
      && check_lower_bound("GradSolverConvergeCheckFrequency", grad_solver_convergence_check_frequency, 0, false)
      && check_nonempty_index_sequence("SourceVertices", source_vertices);
}

template<typename T>
void print_value(const std::string &name, T value)
{
  std::cout << name << ": " << value << std::endl;
}

void Parameters::output_options() {
  std::cout << "======== Solver Parameters ========" << std::endl;
  print_value("HeatSolverMaxIter", heat_solver_max_iter);
  print_value("HeatSolverEps", heat_solver_eps);
  print_value("GradSolverMaxIter", grad_solver_max_iter);
  print_value("GradSolverEps", grad_solver_eps);
  print_value("Penalty", penalty);

  std::cout << "Source vertices: ";
  for(int i = 0; i < static_cast<int>(this->source_vertices.size()); ++ i){
    std::cout << source_vertices[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "====================================" << std::endl;
}
