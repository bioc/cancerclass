\name{cancerclass-internal}
\alias{cancerclass-internal}
\alias{ilogit}
\alias{calc.roc}
\alias{calc.auc}
\alias{get.d}
\alias{get.d2}
\alias{get.prop}
\alias{get.ntrain}
\alias{prepare}
\alias{get.lm}
\alias{filter}
 

\title{ Internal Functions in the cancerclass Package }

\description{
  Internal functions in the cancerclass package,
  which are only user-visible because of the special nature of the cancerclass name space.
}

\usage{
get.d(X, y, method="cor")
get.d2(X, x1, x2, method="cor") 
ilogit(x)
calc.roc(x, y, ci="wilson")
calc.auc(S)
get.prop(n, m, ci="exact", conf.level=.9)
get.ntrain(n1, n2, percent, n.min=5)
prepare(eset, class="class")
get.lm(v, ngenes=10, dist="cor")
filter(method)
 


#.First.lib(lib, pkg)
.onLoad(lib, pkg)
.initFoo(where)
#ratio(x,c1,c2)
}

\arguments{
  \item{X}{ A numeric data matrix or data frame }
  \item{method}{ A numeric data matrix or data frame }
  \item{x}{ Numerical or character vector describing class membership }
  \item{y}{ Fixed number of features }
}

\details{
Internal function are used for auxiliary calculation. 
}

\value{
Each internal function represents function specific results.
}

\author{ Cancerclass Team }
\keyword{internal}
