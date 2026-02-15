# Formatting
- When writing Julia code, try to keep the maximum line length under _92 characters_.
- When writing commit messages, follow the format `component: Brief summary` for
  the title. In the body of the commit message, provide a brief prose summary of
  the purpose of the changes made.
  When referencing external GitHub PRs or issues, use proper GitHub interlinking
  format (e.g., `owner/repo#123` for PRs/issues).

# Coding rules
- When writing functions, avoid generic types of Any if possible.

- For function calls with keyword arguments, use an explicit `;` for clarity.

- For AI agents: **ONLY INCLUDE COMMENTS WHERE TRULY NECESSARY**.
  When the function name or implementation clearly indicates its purpose or
  behavior, redundant comments are unnecessary.

- On the other hand, for general utilities that expected to be used in multiple
  places in this package, it's fine to use docstrings to clarify their
  behavior. However, even in these cases, if the function name and behavior are
  self-explanatory, no special docstring is needed.

# Running test code
Please make sure to test new code that you write.

If explicit test file or code is provided, prioritize running that.
Otherwise, you can run the entire test suite for the TestParticle project by executing
`using Pkg; Pkg.test()` from the root directory of this repository.

# Test code structure
Test code should be written in files that define independent module spaces with
a `test_` prefix (if not already existing in a large `@testset`).
Then include these files from [`test/runtests.jl`](./test/runtests.jl).
This ensures that these files can be run independently from the REPL.
For example, test code for a new feature would be in a file like
this:
> test/test_feature.jl
```julia
module test_feature
using Test # Each module space needs to explicitly declare the code needed for execution
...
end # module test_feature
```
And `test/test_feature.jl` is included from `test/runtests.jl` like this:
> test/runtests.jl
```julia
@testset "TestParticle.jl" begin
   ...
   @testset "feature" include("test_feature.jl")
   ...
end
```

In each test file, you are encouraged to use `@testset "testset name"` to
organize our tests cleanly. For code clarity, unless specifically necessary,
avoid using `using`, `import`, and `struct` definitions  inside `@testset`
blocks, and instead place them at the top level.

Also, you are encouraged to use `let`-blocks to ensure that names aren't
unintentionally reused between multiple test cases.

# Environment-related issues
For AI agents: **NEVER MODIFY [Manifest.toml](./Manifest.toml) BY YOURSELF**.
Agents can modify [Project.toml](./Project.toml) if needed.

# About modifications to code you've written
If you, as an AI agent, add or modify code, and the user appears to have made
further manual changes to that code after your response, please respect those
modifications as much as possible.
For example, if the user has deleted a function you wrote, do not reintroduce
that function in subsequent code generation.
If you believe that changes made by the user are potentially problematic,
please clearly explain your concerns and ask the user for clarification.
