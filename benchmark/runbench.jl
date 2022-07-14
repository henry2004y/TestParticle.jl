# run this script will generate a benchmark report named report.md

using PkgBenchmark

current = BenchmarkConfig(juliacmd = `julia -O3`)
result = benchmarkpkg("TestParticle", current; retune=true)
export_markdown("report.md", result)