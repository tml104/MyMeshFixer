#include <iostream>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <spdlog/spdlog.h>
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "argparser.hpp"
#include "Timer.hpp"

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

namespace MeshUtils
{

	double CalculateAngle(const MyMesh::Point &a, const MyMesh::Point &b, const MyMesh::Point &c)
	{
		MyMesh::Point ab = b - a;
		MyMesh::Point ac = c - a;
		double dot_product = OpenMesh::dot(ab, ac);
		double norm_ab = ab.norm();
		double norm_ac = ac.norm();
		return acos(dot_product / (norm_ab * norm_ac)) * 180.0 / M_PI;
	}

	double CalculateFaceMinAngle(MyMesh &mesh, MyMesh::FaceHandle fh)
	{
		std::vector<MyMesh::Point> points;
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it)
		{
			points.push_back(mesh.point(*fv_it));
		}

		if (points.size() != 3)
		{
			std::cerr << "Error: Face is not a triangle!" << std::endl;
			return 0.0;
		}

		double angles[3];
		for (int i = 0; i < 3; ++i)
		{
			MyMesh::Point v0 = points[i];
			MyMesh::Point v1 = points[(i + 1) % 3];
			MyMesh::Point v2 = points[(i + 2) % 3];

			MyMesh::Point vec1 = v1 - v0;
			MyMesh::Point vec2 = v2 - v0;

			double dot_product = OpenMesh::dot(vec1, vec2);
			double length1 = vec1.norm();
			double length2 = vec2.norm();

			angles[i] = std::acos(dot_product / (length1 * length2));
		}

		double min_angle = std::min({angles[0], angles[1], angles[2]});
		return min_angle * 180.0 / M_PI; // 转换为度数
	}

	double CalculateLength(MyMesh &mesh, MyMesh::HalfedgeHandle heh)
	{
		auto v1 = mesh.from_vertex_handle(heh);
		auto v2 = mesh.to_vertex_handle(heh);
		MyMesh::Point p1 = mesh.point(v1);
		MyMesh::Point p2 = mesh.point(v2);
		MyMesh::Point vec1 = p2 - p1;
		double length = vec1.norm();

		return length;
	}

} // MeshUtils

struct INPUT_ARGUMENTS
{
	std::string file_path;
	double distance_threshold;
	double angle_threshold;
	double scale_factor;
};

struct MeshFixer
{

	double VERTICES_MERGE_DISTANCE_EPSILON = 0.001; // 目前硬编码一个数值
	double FLIP_ANGLE_EPSILON = 150.0;				// 单位：度
	double SCALE_FACTOR = 1.0;

	MyMesh &mesh;

	int flip_count = 0;
	int collapse_count = 0;

	/*
		倍率缩放
	*/
	void ScaleOperation()
	{
		spdlog::info("ScaleOperation start.");
		for (auto vertex_it = mesh.vertices_begin(); vertex_it != mesh.vertices_end(); ++vertex_it)
		{
			MyMesh::Point point = mesh.point(*vertex_it);
			point *= SCALE_FACTOR;
			mesh.set_point(*vertex_it, point);
		}
		spdlog::info("ScaleOperation end.");
	}

	/*
		折叠小于阈值的点
	*/
	void CollapseVerticesOperation()
	{
		spdlog::info("CollapseVerticesOperation start.");

		// 注意：如果要collapse，这几个函数必须调用
		mesh.request_vertex_status();
		mesh.request_edge_status();
		mesh.request_face_status();

		for (MyMesh::HalfedgeIter e_it = mesh.halfedges_begin(); e_it != mesh.halfedges_end(); ++e_it)
		{
			MyMesh::HalfedgeHandle heh = *e_it;

			if (mesh.is_collapse_ok(heh))
			{
				double length = MeshUtils::CalculateLength(mesh, heh);

				if (length < VERTICES_MERGE_DISTANCE_EPSILON) // 小于阈值就折叠
				{
					mesh.collapse(*e_it);
					collapse_count++;
				}
			}
		}

		spdlog::info("CollapseVerticesOperation end.");
	}

	/*
		Flip对角过大的边
		（可能要多次迭代）
	*/
	void FlipEdgesOperation()
	{
		spdlog::info("FlipEdgesOperation start.");

		for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
		{
			try
			{
				MyMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(*e_it, 0);
				MyMesh::HalfedgeHandle heh1 = mesh.halfedge_handle(*e_it, 1); // 等等，是非流形怎么办？

				MyMesh::VertexHandle v0 = mesh.to_vertex_handle(heh0);
				MyMesh::VertexHandle v1 = mesh.to_vertex_handle(heh1);

				MyMesh::VertexHandle v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh0));
				MyMesh::VertexHandle v3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh1));

				double angle0 = MeshUtils::CalculateAngle(mesh.point(v2), mesh.point(v1), mesh.point(v0));
				double angle1 = MeshUtils::CalculateAngle(mesh.point(v3), mesh.point(v1), mesh.point(v0));

				if (angle0 > FLIP_ANGLE_EPSILON || angle1 > FLIP_ANGLE_EPSILON)
				{
					spdlog::info("angle0: {}, angle1: {}", angle0, angle1);

					if (mesh.is_flip_ok(*e_it))
					{
						mesh.flip(*e_it);
						spdlog::info("fliped. flip count: {}", flip_count);
						flip_count++;
					}
				}
			}
			catch (const std::exception &e)
			{
				spdlog::error("Caught an exception in FlipEdgesOperation: {}", e.what());
			}
			catch (...)
			{
				spdlog::error("Caught an unknown exception in FlipEdgesOperation");
			}
		}

		spdlog::info("FlipEdgesOperation end.");
	}

	/*
		以直方图形式打印边长统计信息
	*/
	void PrintStatsForEdgeLength()
	{
		spdlog::info("PrintStatsForEdgeLength start.");

		std::vector<double> lengths;

		for (MyMesh::HalfedgeIter e_it = mesh.halfedges_begin(); e_it != mesh.halfedges_end(); ++e_it)
		{
			MyMesh::HalfedgeHandle heh = *e_it;

			double length = MeshUtils::CalculateLength(mesh, heh);
			lengths.emplace_back(length);
		}

		std::sort(lengths.begin(), lengths.end());

		double min_length = lengths.front();
		double max_length = lengths.back();

		int part = 10;
		
		double delta = (max_length - min_length) / part;
		double current_l = min_length;
		double current_r = min_length + delta;
		int cnt = 0;
		
		auto print_line = [&](){
			spdlog::info("length_range [{} - {}): {} ({})", current_l, current_r, cnt, cnt*1.0/lengths.size());
		};

		for(auto l: lengths){
			if(l<current_r){
				cnt++;
			}else{
				print_line();
				current_l +=delta;
				current_r +=delta;
				cnt = 0;
			}
		}
		print_line();
		spdlog::info("min_length: {}", min_length);
		spdlog::info("max_length: {}", max_length);

		spdlog::info("PrintStatsForEdgeLength end.");
	}

	void PrintMeshInfo()
	{
		int n_vertices = mesh.n_vertices();
		int n_edges = mesh.n_edges();
		int n_faces = mesh.n_faces();

		spdlog::info("=== Mesh Info ===");
		spdlog::info("n_vertices: {}", n_vertices);
		spdlog::info("n_edges: {}", n_edges);
		spdlog::info("n_faces: {}", n_faces);
		spdlog::info("=== Mesh Info End ===");
	}

	void Status()
	{
		spdlog::info("=== Status for Mesh Operation ===");
		spdlog::info("collapse count: {}", collapse_count);
		spdlog::info("flip count: {}", flip_count);
		spdlog::info("=== Status for Mesh Operation End ===");
	}

	void Start(INPUT_ARGUMENTS input_args)
	{
		spdlog::info("=== Before Operations ===");

		VERTICES_MERGE_DISTANCE_EPSILON = input_args.distance_threshold;
		FLIP_ANGLE_EPSILON = input_args.angle_threshold;
		SCALE_FACTOR = input_args.scale_factor;

		PrintMeshInfo();
		PrintStatsForEdgeLength();

		if (SCALE_FACTOR != 1.0)
		{
			ScaleOperation();
		}

		CollapseVerticesOperation();
		mesh.garbage_collection();
		FlipEdgesOperation();

		spdlog::info("=== After Operations ===");
		PrintMeshInfo();
		Status();
	}

	void Clear()
	{
		flip_count = 0;
		collapse_count = 0;
	}

	MeshFixer(MyMesh &mesh) : mesh(mesh) {}
};

void Start(INPUT_ARGUMENTS input_args)
{
	spdlog::set_pattern("[%H:%M:%S %z] [%n] [%^%L%$] [thread %t] [%s] [%@] %v");

	std::string output_file_path = input_args.file_path.substr(0, input_args.file_path.length() - 4) + "_output.obj";

	//std::string file_name = "B_ent2(1).stl";
	//std::string output_file_name = "B_ent2(1)_output.stl";

	Timer::Timer timer;

	MyMesh mesh;

	if (!OpenMesh::IO::read_mesh(mesh, input_args.file_path))
	{
		spdlog::error("Cannot read mesh from: {}", input_args.file_path);
		return;
	}
	timer.Split("Input Mesh");

	MeshFixer meshFixer(mesh);
	meshFixer.Start(input_args);
	timer.Split("Process");

	if (!OpenMesh::IO::write_mesh(mesh, output_file_path))
	{
		spdlog::error("Cannot write mesh from: {}", output_file_path);
		return;
	}
	timer.Split("Output Mesh");

	timer.PrintTimes();
}

int main(int argc, char const *argv[])
{
	// parse args
	auto args_parser = util::argparser("ObjFixer by TML104");
	args_parser.set_program_name("MyMeshFixer")
		.add_help_option()
		.use_color_error()
		.add_argument<std::string>("input_obj_file", "stl model path")
		.add_option<double>("-d", "--distance", "distance threshold", 0.001)
		.add_option<double>("-a", "--angle", "angle threshold", 150)
		.add_option<double>("-s", "--scale", "scale factor", 1.0)
		.parse(argc, argv);

	INPUT_ARGUMENTS input_args;
	input_args.distance_threshold = args_parser.get_option_double("-d");
	input_args.angle_threshold = args_parser.get_option_double("-a");
	input_args.scale_factor = args_parser.get_option_double("-s");
	input_args.file_path = args_parser.get_argument<std::string>("input_obj_file");

	Start(input_args);
	return 0;
}