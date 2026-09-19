// E9 Phase-1 unit test: ExtTriangleMesh vertex-motion series
#include <cstdio>
#include <cmath>
#include <vector>
#include "luxrays/core/exttrianglemesh.h"
#include "luxrays/core/geometry/transform.h"

using namespace luxrays;

static int fails = 0;
#define CHECK(cond, msg) do { if(!(cond)) { printf("FAIL: %s\n", msg); ++fails; } else printf("PASS: %s\n", msg); } while(0)

static ExtTriangleMeshUPtr makeQuad(float z0, float z1) {
	VertexBuffer vs(4);
	vs[0] = Point(0,0,0); vs[1] = Point(1,0,0);
	vs[2] = Point(1,1,0); vs[3] = Point(0,1,0);
	TriangleBuffer ts(2);
	ts[0] = Triangle(0,1,2); ts[1] = Triangle(0,2,3);
	NormalBuffer ns;  // no normals

	auto m = std::make_unique<ExtTriangleMesh>(std::move(vs), std::move(ts), std::move(ns));

	std::vector<VertexBuffer> steps;
	steps.emplace_back(4); steps.emplace_back(4);
	for (u_int v = 0; v < 4; ++v) { steps[0][v] = Point(v%2, v/2, z0); }
	for (u_int v = 0; v < 4; ++v) { steps[1][v] = Point(v%2, v/2, z1); }
	m->SetVertexMotion(std::vector<float>{0.f, 1.f}, std::move(steps));
	return m;
}

int main() {
	auto m = makeQuad(0.f, 2.f);

	// basics
	CHECK(m->HasVertexMotion(), "HasVertexMotion true");
	CHECK(m->GetVertexMotionStepCount() == 2, "2 steps");
	CHECK(m->GetVertexMotionTimes().size() == 2, "2 times");

	// clamping + lerp
	Point p0 = m->GetVertexAtTime(0, -1.f);
	CHECK(fabsf(p0.z - 0.f) < 1e-6, "clamp below -> step0");
	Point p1 = m->GetVertexAtTime(0, 2.f);
	CHECK(fabsf(p1.z - 2.f) < 1e-6, "clamp above -> step1");
	Point pm = m->GetVertexAtTime(0, 0.25f);
	CHECK(fabsf(pm.z - 0.5f) < 1e-6, "lerp 0.25 -> z=0.5");

	// no motion -> static vertex
	{
		VertexBuffer vs2(1); vs2[0] = Point(9,9,9);
		TriangleBuffer ts2; // empty
		// build a valid 1-tri mesh instead
		VertexBuffer v3(3);
		v3[0]=Point(0,0,0); v3[1]=Point(1,0,0); v3[2]=Point(0,1,0);
		TriangleBuffer t3(1); t3[0]=Triangle(0,1,2);
		NormalBuffer n3;
		ExtTriangleMesh plain(std::move(v3), std::move(t3), std::move(n3));
		CHECK(!plain.HasVertexMotion(), "plain mesh no motion");
		Point ps = plain.GetVertexAtTime(1, 0.5f);
		CHECK(fabsf(ps.x-1.f)<1e-6, "no-motion -> static vertex");
	}

	// validation errors
	{
		VertexBuffer v3(3); v3[0]=Point(0,0,0); v3[1]=Point(1,0,0); v3[2]=Point(0,1,0);
		TriangleBuffer t3(1); t3[0]=Triangle(0,1,2);
		NormalBuffer n3;
		ExtTriangleMesh bad(std::move(v3), std::move(t3), std::move(n3));
		bool threw = false;
		try {
			std::vector<VertexBuffer> s2; s2.emplace_back(3); s2.emplace_back(2); // wrong count
			bad.SetVertexMotion(std::vector<float>{0.f,1.f}, std::move(s2));
		} catch (...) { threw = true; }
		CHECK(threw, "vert-count mismatch throws");
		threw = false;
		try {
			std::vector<VertexBuffer> s2; s2.emplace_back(3); s2.emplace_back(3);
			bad.SetVertexMotion(std::vector<float>{1.f,0.f}, std::move(s2)); // non-monotonic
		} catch (...) { threw = true; }
		CHECK(threw, "non-monotonic throws");
		threw = false;
		try {
			std::vector<VertexBuffer> s2; s2.emplace_back(3);
			bad.SetVertexMotion(std::vector<float>{0.f}, std::move(s2)); // <2 steps
		} catch (...) { threw = true; }
		CHECK(threw, "single step throws");
	}

	// CopyExt preserves the series (no vertex override)
	{
		auto c = m->CopyExt(std::nullopt, std::nullopt, std::nullopt,
				std::nullopt, std::nullopt, std::nullopt, 0.f);
		CHECK(c->HasVertexMotion(), "CopyExt keeps motion");
		CHECK(fabsf(c->GetVertexAtTime(0, 0.5f).z - 1.f) < 1e-6, "copy lerp z=1");
	}
	// CopyExt with overridden vertices drops the series
	{
		VertexBuffer nv(4);
		for (u_int i = 0; i < 4; ++i) nv[i] = Point(i,i,i);
		auto c = m->CopyExt(std::move(nv), std::nullopt, std::nullopt,
				std::nullopt, std::nullopt, std::nullopt, 0.f);
		CHECK(!c->HasVertexMotion(), "CopyExt drops motion on vertex override");
	}

	// ApplyTransform transforms step buffers
	{
		auto c = m->Copy();
		Matrix4x4 mat = Matrix4x4::MAT_IDENTITY;
		mat.m[0][3] = 10.f; // translate +10 x
		c->ApplyTransform(Transform(mat));
		Point p = c->GetVertexAtTime(0, 1.f);
		CHECK(fabsf(p.x - 10.f) < 1e-6 && fabsf(p.z - 2.f) < 1e-6,
				"ApplyTransform moves step verts");
	}

	// Merge: identical times -> merged series; partial presence -> throw
	{
		auto a = makeQuad(0.f, 2.f);
		auto b = makeQuad(5.f, 7.f);
		std::vector<std::reference_wrapper<const ExtTriangleMesh>> in = {*a, *b};
		auto merged = ExtTriangleMesh::Merge(in, std::nullopt);
		CHECK(merged->HasVertexMotion(), "Merge keeps motion (identical times)");
		CHECK(merged->GetTotalVertexCount() == 8, "merged 8 verts");
		CHECK(fabsf(merged->GetVertexAtTime(5, 1.f).z - 7.f) < 1e-6,
				"merged step offset correct");

		// partial presence -> throw
		VertexBuffer v3(3); v3[0]=Point(0,0,0); v3[1]=Point(1,0,0); v3[2]=Point(0,1,0);
		TriangleBuffer t3(1); t3[0]=Triangle(0,1,2);
		NormalBuffer n3;
		auto plain = std::make_unique<ExtTriangleMesh>(std::move(v3), std::move(t3), std::move(n3));
		std::vector<std::reference_wrapper<const ExtTriangleMesh>> in2 = {*a, *plain};
		bool threw = false;
		try { ExtTriangleMesh::Merge(in2, std::nullopt); } catch (...) { threw = true; }
		CHECK(threw, "Merge throws on partial motion");

		// different times -> throw
		auto c = makeQuad(0.f, 2.f);
		// rebuild c with different times
		{
			std::vector<VertexBuffer> st;
			st.emplace_back(4); st.emplace_back(4);
			for (u_int v = 0; v < 4; ++v) { st[0][v]=Point(0,0,0); st[1][v]=Point(0,0,1); }
			c->SetVertexMotion(std::vector<float>{0.f, 0.7f}, std::move(st));
		}
		std::vector<std::reference_wrapper<const ExtTriangleMesh>> in3 = {*a, *c};
		threw = false;
		try { ExtTriangleMesh::Merge(in3, std::nullopt); } catch (...) { threw = true; }
		CHECK(threw, "Merge throws on different times");
	}

	// serialization drops motion (static fallback)
	{
		m->SaveSerialized("/tmp/e9_test/motion.bpy");
		auto loaded = ExtTriangleMesh::LoadSerialized("/tmp/e9_test/motion.bpy");
		CHECK(!loaded->HasVertexMotion(), "serialized mesh drops motion (static fallback)");
		CHECK(loaded->GetTotalVertexCount() == 4, "serialized verts intact");
	}

	printf("\n%d failures\n", fails);
	return fails ? 1 : 0;
}
