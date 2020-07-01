<!-- =============================
     ABOUT
    ============================== -->

\begin{:section, title="About this Package", name="About"}

![CompatHelper](https://github.com/fusion809/FunctionIntegrator.jl/workflows/CompatHelper/badge.svg)
![TagBot](https://github.com/fusion809/FunctionIntegrator.jl/workflows/TagBot/badge.svg)
![Travis](https://travis-ci.com/fusion809/FunctionIntegrator.jl.svg?branch=master)

[FunctionIntegrator.jl](https://github.com/fusion809/FunctionIntegrator.jl) should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not. Likewise this package is less user-friendly as it requires you to decide which $N$ value you are going to go with.

This package provides the following functions:

* `chebyshev_quadrature(f::Function, N::Number, k::Integer, a::Number, b::Number)`
* `hermite_quadrature(f::Function, N::Number, k::Integer)`
* `jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)`
* `laguerre_quadrature(f::Function, N::Number, k::Integer)`
* `legendre_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `lobatto_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `radau_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule_left(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule_midpoint(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule_right(f::Function, N::Number, a::Number, b::Number)`
* `simpsons_rule(f::Function, N::Number, a::Number, b::Number)`
* `trapezoidal_rule(f::Function, N::Number, a::Number, b::Number)`

use Julia's help function (e.g. by typing `?chebyshev_quadrature`) to find out usage information, should you need it. The choice of function table below also explains some of the details of each of these functions, such as their arguments.

\end{:section}

<!-- =============================
     GETTING STARTED
     ============================== -->
\begin{:section, title="Getting started"}
This package is currently in Julia's General registry, thus to install it one merely needs to run:

```julia-repl
(v1.4) pkg> add FunctionIntegrator
```

and import it using:

```julia-repl
julia> using FunctionIntegrator
```
\end{:section}

<!-- =============================
     SPECIAL COMMANDS
     ============================== -->
\begin{:section, title="Choice of function"}
As a general rule of thumb, `simpsons_rule` should be the function you use when you are unsure which function to use as its approximations with large N (e.g. $1 \times 10^{4}$) are nearly always accurate to at least 6 digits. Despite this, for many problems some of the `_quadrature` functions may provide more accurate results with a smaller N value. The main time when `simpsons_rule` should be avoided is when there are unremovable singularities at the endpoints of the domain of integration, in which case using `chebyshev_quadrature` with $k=1$ or using `legendre_quadrature` is likely best.

Each of the functions whose name ends with `_quadrature` uses [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature), the specifics of which differ between functions. Out of them, `legendre_quadrature` is perhaps the best to go with when you are uncertain which out of the `_quadrature` functions to go with as it can be applied to any grid and any problem and *usually* arrives at a result that is accurate to at least 7 decimal places provided $N \geq 1 \times 10^{4}$. `legendre_quadrature` has the added benefit of being perhaps the fastest `_quadrature` function after controlling for accuracy of the result obtained.

The [test/](https://github.com/fusion809/FunctionIntegrator.jl/tree/master/test/) folder has test scripts that approximate various different integrals (each file, except [test/runtests.jl](https://github.com/fusion809/FunctionIntegrator.jl/blob/master/test/runtests.jl), pertains to a different integral); the N values given are the smallest possible to pass each of the tests listed (except when the test involves a less than (<) sign). If you want to know which function to use for which integral, these tests may be useful as a rough guide.

~~~
<figure float="center">
    <img src="/assets/Root_mean_square_computation_time_FunctionIntegrator.jl.svg" width="100%">
    <figcaption><b>Figure 1: a bar graph to compare the computation times for each fo the flexible-domain methods provided by the package, except <code>rectange_rule_left</code> and <code>rectangle_rule_right</code>.<sup>1</sup> Data is from <a href="https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/355708079" link="_blank">this build</a>.</b></caption>
</figure>
~~~

| Function               | Domain~~~<sup>2</sup>~~~ | Weight | Arguments | Notes |
|------------------------|--------|--------|-----------|-------|
| `chebyshev_quadrature` | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | $k=1$: $\dfrac{1}{\sqrt{1-x^2}}$~~~<br/><br/>~~~$k=2$: $\sqrt{1-x^2}$~~~<br/><br/>~~~$k=3$: $\sqrt{\dfrac{1+x}{1-x}}$~~~<br/><br/>~~~$k=4$: $\sqrt{\dfrac{1-x}{1+x}}$ | `f`, the function being integrated. ~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k`, the type of Chebyshev quadrature being used. 1, 2, 3, and 4 refer to the Chebyshev $T_n$, $U_n$, $V_n$ and $W_n$ polynomials respectively.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Chebyshev-Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) (note this article does not mention 3rd and 4th quadrature types, corresponding to $k=3$ and $k=4$, respectively). If there are unremovable singularities at the endpoints, it with $k=1$ or `legendre_quadrature` are preferred. |
| `hermite_quadrature`   | $[-\infty,\infty]$ | $e^{-x^2}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x^2}$ is assumed to be part of the integrand ($k=2$) or not). | Uses [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature). Only use this if your integration domain is $[-\infty,\infty]$ and your integrand rapidly goes to zero as the absolute value of $x$ gets larger. |
| `jacobi_quadrature`    | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | $(1-x)^{\alpha}(1+x)^{\beta}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`α` and `β` are parameters of the weighting function.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature). After controlling for `N` it is one of the slowest quadrature methods. A condition of Gauss-Jacobi quadrature is that $\alpha, \beta \gt -1$. When $\alpha = \beta = 0$, this reduces to `legendre_quadrature`. |
| `laguerre_quadrature`  | $[0,\infty]$ | $e^{-x}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x}$ is assumed to be part of the integrand ($k=2$) or not). | Uses [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature). Only use this if your integration domain is $[0,\infty]$ and your integrand rapidly goes to $0$ as $x$ gets larger. |
| `legendre_quadrature`  | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature). Generally, this is the best `_quadrature` function to go with when you are otherwise unsure which to go with. |
| `lobatto_quadrature`   | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Gauss-Lobatto quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Lobatto_rules). This function includes, in the calculation, the values of the integrand at one of the endpoints. Consequently, if there are unremovable singularities at the endpoints, this function may fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `radau_quadrature`     | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses Gauss–Radau quadrature, for which there is no Wikipedia article, the best article (simplest) I could find on it are [these lecture notes](https://web.archive.org/web/20200628202423/http://www.math.hkbu.edu.hk/ICM/LecturesAndSeminars/08OctMaterials/2/Slide3.pdf). This function includes, in the calculation, the values of the function at the endpoints. Consequently, if there are unremovable singularities at either or both of the endpoints, this function will fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `rectangle_rule_left`  | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses the rectangle rule, specifically the [left Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum#Left_Riemann_sum). Usually this or `rectangle_rule_right` is the least accurate method. In fact, many of the tests in the FunctionIntegrator.jl repository fail to get accuracy to 7 significant figures with `rectangle_rule_left` with any practically viable value of `N`. |
| `rectangle_rule_midpoint`  | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses the rectangle rule, specifically the [Riemann midpoint rule](https://en.wikipedia.org/wiki/Riemann_sum#Midpoint_rule). Usually this is more accurate than `rectangle_rule_left` and `rectangle_rule_right` and sometimes rivals `trapezoidal_rule` for accuracy. Interestingly, going by my Travis tests it appears to be even more efficient than `simpsons_rule`. |
| `rectangle_rule_right`  | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses the rectangle rule, specifically the [right Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum#Right_Riemann_sum). Usually this or `rectangle_rule_left` is the least accurate method. In fact, many of the tests in the FunctionIntegrator.jl repository fail to get accuracy to 7 significant figures with `rectangle_rule_right` with any practically viable value of `N`. |
| `simpsons_rule`        | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule). It is one of the best functions to use when you are unsure which to use, provided there are no unremovable singularities within the integration domain, including the endpoints. |
| `trapezoidal_rule`     | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule). It has the same caveats as `simpsons_rule`. |

**Notes**:
1. `rectangle_rule_left`, `rectangle_rule_midpoint` and `rectangle_rule_right` are not included because they failed to provide the required level of accuracy for many tests with all practically viable `N` values.
2. The $x\in[-1,1]$ refers to the quadrature nodes, which are also referred to in the weighting function column.
\end{:section}


<!-- =============================
     SHOWING CODE
    ============================== -->

\begin{:section, title="Showing Code"}

\lead{
    Franklin can run your Julia code on the fly and show the output.
}

**Setting up the environment**: the first step is to ensure that the folder with your source has the proper environment including your package.
To do so, in the Julia REPL, navigate to the source (e.g. `cd("page/")`), activate the environment (e.g. `using Pkg; Pkg.activate()`) and add the package(s) that you need (e.g. `Pkg.add("DataFrames")`).
If you check the status or the Project.toml, you will see that `Franklin` is already in there on top of whatever packages you might have chosen to add.
In our current case:

```
Status `~/.julia/dev/PkgPage/page/Project.toml`
  [a93c6f00] DataFrames v0.21.2
  [713c75ef] Franklin v0.8.2
```

Also add the package in the `DeployPage.yml` e.g. in our case there is:

```julia
Pkg.add(["NodeJS", "DataFrames"]);
```

Once that's set up, you can use "named" code blocks i.e. code blocks that look like

`````
```julia:ex
using DataFrames
df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
first(df, 3)
```
`````

where `:ex` is the "named part" (`ex` being the name, which should be unique on the page).

```julia:ex
using DataFrames
df = DataFrame(A = 1:4, B = ["M", "F", "F", "M"])
first(df, 3)
```

You can control the indentation and appearance of the output block in the `config.md` too.

\end{:section}


<!-- =============================
     Deploying
    ============================== -->

\begin{:section, title="Deployment"}

\lead{Make your page available online easily by leveraging GitHub Actions and GitHub Pages.}

By following these instructions, the content of the rendered website will be copied to a `gh-pages` branch where it will be deployed by GitHub.
If you would like to deploy the page with your own URL or using something else than GitHub, have a look at the specific instructions further on.

**Adjust DeployPage**: start by checking the `.github/workflows/DeployPage.yml` in particular:
* if you want to use Python or matplotlib, uncomment the relevant lines
* in the `run` block ensure that
    * `NodeJS` and `PkgPage` are added,
    * any packages that your page might rely on are added,
    * the `optimize` call has the appropriate `input` and `output` path (if you're in the default setting, leave as is).

**Keys**: in order to have your page be built and deployed on GitHub, you will need to generate a keypair and add it to the GitHub repo. To do so:

1. run in your terminal `ssh-keygen -N "" -f franklin`,
1. copy the entire content of the generated `franklin` file and put it as a new secret named `FRANKLIN_PRIV` on <https://github.com/USERNAME/PACKAGE.jl/settings/secrets/new>,
1. copy the entire content of the generated `franklin.pub` file and put it as a new deploy key named `FRANKLIN_PUB` on <https://github.com/USERNAME/PACKAGE.jl/settings/keys> with `read/write` access,
1. remove both files.

**GitIgnore**: it's important you specify that `page/__site` should be ignored by git and not pushed to your repository otherwise the build process might not work properly. To do so create a file `.gitignore` containing the line

```
page/__site
```

as shown [here](https://github.com/tlienart/PkgPage.jl/blob/cce098535eb95c2c3ba919d605792abfee57710c/.gitignore#L3).

**GitAttributes**: in order for GitHub to ignore `page` folder it the language statistics for your repository, make sure to add a file `.gitattributes` with content

```
page/* linguist-vendored
```

like [this](https://github.com/tlienart/PkgPage.jl/blob/master/.gitattributes).

Now if you push your changes and, generally, whenever the `master` branch of your package gets updated, the  build process will be triggered and your page updated and deployed.
**That's it**.

**Avoiding clashes with Documenter.jl**: if you already use [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) you might want your page to be deployed in a specific folder of `gh-pages` as Documenter also generates files in `gh-pages`.

\alert{This will typically not be necessary as the names created by PkgPage and Documenter don't clash, but you might still prefer to avoid mixing the two (in which case, read on).}

you can do so in two steps:

1. change the `run` part of `DeployPage.yml` by specifying the `output` keyword argument  in `PkgPage.optimize` for instance: `PkgPage.optimize(input="page", output="page")`,
1. change the `prepath` in `config.md` to reflect that the base URL will contain that additional folder, for instance `@def prepath = "YourPackage.jl/page"`.

**Use your own URL**: you can usually get host services like Netlify to deploy a specific branch of a GitHub repo, do make sure to set `@def prepath = ""` in your `config.md` though.

If you want to do the deployment without GitHub actions then you will need to:

* ensure you have `purgecss` and `highlights` installed and available to `NodeJS`, the simplest way to do this is to install them via `NodeJS` with

```
using NodeJS;
run(`$(npm_cmd()) install highlight.js`);
run(`$(npm_cmd()) install purgecss`);
```
\\
* run `PkgPage.optimize(input="page", output="")` (adapting `input` as required)
* place the content of `page/__site` wherever your server requires it.

\end{:section}
