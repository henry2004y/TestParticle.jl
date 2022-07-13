using PkgBenchmark

current = BenchmarkConfig(juliacmd = `julia -O3`)
result = benchmarkpkg("TestParticle", current)
export_markdown("report.md", result)