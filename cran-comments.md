## R CMD check results

0 errors | 0 warnings | 1 note

This is a new release, being resubmitted following an initial review.

We are grateful to CRAN for the manual review of the package and feedback provided. We have made changes as follows:

> Warning: Examples in comments in:
>   delete_built_in_proposal.Rd
>
> The reason for commenting out is clear, however, wrapping it in
> \dontrun{} is the better solution here if no one should run it. If they
> just should not be tested, \donttest{} would be ideal.

We have uncommented the examples in `delete_built_in_proposal.Rd` and wrapped them in `\dontrun{}` as requested.

> Please ensure that your functions do not write by default or in your
> examples/vignettes/tests in the user's home filespace (including the
> package directory and getwd()). This is not allowed by CRAN policies.
> Please omit any default path in writing functions. In your
> examples/vignettes/tests you can write to tempdir(). -> R/zzz.R

We have removed the automatic creation of grid files which used to take place during `.onLoad` in `R/zzz.R`. Now it will only attempt to read any grid files that already exist. We have therefore changed to use pre-generated default grid files which we distribute in the package `inst/` folder, to avoid the need to generate them on first package load.

To the best of our understanding of the review, we believe this fully addresses the request: please let us know if the intention was to cover anything else. For context which may help the review: the only other writing we do is to the directory obtained from `tools::R_user_dir()`, and those remaining cases are only when the user explicitly runs a function which is clearly named and documented to trigger writing of configuration/grid caches to that directory. We believe this to be in line with the CRAN Repository Policy which states:

"For R version 4.0 or later (hence a version dependency is required or only conditional use is possible), packages may store user-specific data, configuration and cache files in their respective user directories obtained from tools::R_user_dir(), provided that by default sizes are kept as small as possible and the contents are actively managed (including removing outdated material)."

With many thanks for the time and effort in reviewing our package.
 