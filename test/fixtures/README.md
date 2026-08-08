# Golden Fixtures

`generate_golden.R` evaluates the local R source with deterministic point sets
and writes CSV fixtures for the public geometry pipeline. Regenerate them with:

```powershell
Rscript test/fixtures/generate_golden.R R test/fixtures/generated
```

The generated data is the compatibility oracle for the Julia implementation.
