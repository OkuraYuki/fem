#include "electrostatic_analyzer.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#endif

struct ModelSpec {
	std::string name;
	std::string mesh_file;
	std::string material_file;
	std::string initial_file;
	std::string pulse_file;
	std::string config_file;
	std::string output_dir;
};

static bool file_exists(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) return false;
	fclose(fp);
	return true;
}

static std::string trim(const std::string &text)
{
	size_t start = text.find_first_not_of(" \t\r\n");
	if (start == std::string::npos) return std::string();
	size_t end = text.find_last_not_of(" \t\r\n");
	return text.substr(start, end - start + 1);
}

static std::string join_path(const std::string &lhs, const std::string &rhs);  // Forward declaration

static ModelSpec make_model_spec(const std::string &name)
{
	ModelSpec spec;
	spec.name = name;
	spec.mesh_file = join_path("models", join_path(spec.name, "in.dat"));
	spec.material_file = join_path("models", join_path(spec.name, "sin.dat"));
	spec.initial_file = join_path("models", join_path(spec.name, "tmate.dat"));
	spec.pulse_file = join_path("models", join_path(spec.name, "pulse.dat"));
	spec.config_file = join_path("models", join_path(spec.name, "sina.dat.t"));
	spec.output_dir = join_path("output", spec.name);
	return spec;
}

static std::vector<ModelSpec> load_model_specs(const char *filename)
{
	std::vector<ModelSpec> specs;
	std::ifstream ifs(filename);
	if (!ifs) return specs;

	std::string line;
	while (std::getline(ifs, line)) {
		line = trim(line);
		if (line.empty() || line[0] == '#') continue;

		ModelSpec spec = make_model_spec(line);

		if (!file_exists(spec.mesh_file.c_str())) {
			printf("Model '%s': mesh file not found at %s\n", spec.name.c_str(), spec.mesh_file.c_str());
			continue;
		}

		specs.push_back(spec);
	}

	return specs;
}

static bool path_exists(const std::string &path)
{
	struct stat st;
	return !path.empty() && stat(path.c_str(), &st) == 0 && (st.st_mode & S_IFDIR);
}

static bool create_directory_recursive(const std::string &path)
{
	if (path.empty() || path_exists(path)) return true;

	size_t pos = path.find_last_of("/\\");
	if (pos != std::string::npos) {
		std::string parent = path.substr(0, pos);
		if (!parent.empty() && !path_exists(parent) && !create_directory_recursive(parent)) return false;
	}

#ifdef _WIN32
	if (_mkdir(path.c_str()) != 0 && !path_exists(path)) return false;
#else
	if (mkdir(path.c_str(), 0755) != 0 && !path_exists(path)) return false;
#endif
	return true;
}

static std::string join_path(const std::string &lhs, const std::string &rhs)
{
	if (lhs.empty()) return rhs;
	if (rhs.empty()) return lhs;
	if (lhs.back() == '/' || lhs.back() == '\\') return lhs + rhs;
	return lhs + "/" + rhs;
}

static int run_one_model(const ModelSpec &spec, bool batch_mode, std::ofstream *batch_csv)
{
	ElectrostaticAnalyzer analyzer;
	std::string output_dir = spec.output_dir.empty() ? spec.name : spec.output_dir;
	if (batch_mode) {
		if (!create_directory_recursive(output_dir)) {
			printf("Cannot create output directory: %s\n", output_dir.c_str());
			return 1;
		}
		const std::string step_prefix = join_path(output_dir, "tempa.dat");
		const std::string report_path = join_path(output_dir, "analysis_summary.txt");
		analyzer.set_output_paths(step_prefix, report_path);
	} else {
		analyzer.set_output_paths("tempa.dat", "analysis_summary.txt");
	}

	printf("=== Model: %s ===\n", spec.name.c_str());
	analyzer.read_mesh(spec.mesh_file.c_str());
	if (analyzer.get_num_nodes() <= 0 || analyzer.get_num_elements() <= 0) {
		printf("Model '%s' was not loaded because the mesh is unavailable or too large.\n", spec.name.c_str());
		return 1;
	}
	analyzer.read_material_and_bc(spec.material_file.c_str());
	analyzer.read_initial_potential(spec.initial_file.c_str());
	if (file_exists(spec.pulse_file.c_str())) {
		analyzer.read_pulse_config(spec.pulse_file.c_str());
	}
	analyzer.read_analysis_config(spec.config_file.c_str());
	analyzer.solve();

	if (batch_csv && batch_csv->is_open()) {
		(*batch_csv)
			<< spec.name << ','
			<< spec.mesh_file << ','
			<< spec.material_file << ','
			<< spec.initial_file << ','
			<< spec.config_file << ','
			<< output_dir << ','
			<< analyzer.get_num_nodes() << ','
			<< analyzer.get_num_elements() << ','
			<< analyzer.get_boundary_count() << ','
			<< analyzer.get_fixed_material_count() << ','
			<< analyzer.get_material_count() << ','
			<< analyzer.get_final_step() << ','
			<< analyzer.get_final_time() << ','
			<< analyzer.get_final_max_diff() << ','
			<< (analyzer.get_converged_early() ? "yes" : "no")
			<< '\n';
	}

	return 0;
}

int main(int argc, char **argv)
{
	const char *manifest = (argc >= 2) ? argv[1] : nullptr;
	const char *model_filter = (argc >= 3) ? argv[2] : nullptr;

	if (manifest) {
		if (argc == 2 && !file_exists(manifest)) {
			ModelSpec spec = make_model_spec(manifest);
			if (!file_exists(spec.mesh_file.c_str())) {
				printf("Cannot open manifest '%s' and model '%s' was not found at %s\n",
					manifest, spec.name.c_str(), spec.mesh_file.c_str());
				return 1;
			}
			return run_one_model(spec, true, nullptr);
		}

		std::vector<ModelSpec> specs = load_model_specs(manifest);
		if (specs.empty()) {
			printf("Cannot open or empty manifest: %s\n", manifest);
			return 1;
		}

		std::ofstream batch_csv("batch_summary.csv");
		if (batch_csv.is_open()) {
			batch_csv << "model_name,mesh_file,material_file,initial_file,config_file,output_dir,nodes,elements,boundaries,fixed_materials,materials,final_step,final_time,final_max_diff,converged_early\n";
		}

		int matched = 0;
		int ran = 0;
		int failed = 0;
		for (size_t i = 0; i < specs.size(); ++i) {
			if (model_filter && specs[i].name != model_filter) continue;
			++matched;
			if (run_one_model(specs[i], true, batch_csv.is_open() ? &batch_csv : nullptr) == 0) {
				++ran;
			} else {
				++failed;
			}
		}

		if (matched == 0) {
			printf("No model matched the filter.\n");
			return 1;
		}

		if (failed > 0) {
			printf("Batch completed with %d failed model(s).\n", failed);
			return 1;
		}

		printf("Wrote batch_summary.csv\n");
		return 0;
	}

	// Fallback: current single-model workflow
	ModelSpec spec;
	spec.name = "default";
	spec.mesh_file = "in.dat";
	spec.material_file = "sin.dat";
	spec.initial_file = "tmate.dat";
	spec.pulse_file = "pulse.dat";
	spec.config_file = file_exists("sina.dat.t") ? "sina.dat.t" : "sina.dat";
	return run_one_model(spec, false, nullptr);
}
