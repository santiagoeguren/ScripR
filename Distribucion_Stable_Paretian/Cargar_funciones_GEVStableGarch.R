

#Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                           DESCRIPTION:
#  'GEVSTABLEGARCH'                  GEVSTABLEGARCH Class representation
################################################################################


# Class Representation:
setClass("GEVSTABLEGARCH",
         representation(
           call = "call",
           formula = "formula",
           method = "character",
           convergence = "numeric",
           messages = "list",
           data = "numeric",
           fit = "list",
           residuals = "numeric",
           h.t = "numeric",
           sigma.t = "numeric",
           title = "character",
           description = "character")
)


################################################################################




# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:               SPECIFICATION:
#  GEVSTABLEGARCHSPEC     S4 GEVSTABLEGARCHSPEC Class representation
################################################################################


# S4 GEVStableGarchSpec Class
setClass("GEVSTABLEGARCHSPEC",
         representation(
           call = "call",
           formula = "formula",
           model = "list",
           presample = "matrix",
           distribution = "character",
           rseed = "numeric")
)


################################################################################


# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:              t3 (now called GAt) distribution proposed in Paolella (1997)
#  pgat                   Probability function for the GAt
#  dgat                   Density for the GAt-distribution
#  qgat                   Quantile function for the GAt
#  rgat                   Random Number Generator for the GAt
################################################################################


dgat <-
  function(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1, log = FALSE)
  {
    # A function implemented by Thiago Sousa

    # Description:
    #   Compute the density for the
    #   so called t3-distribution (now called GAt) defined in Paolella (1997).
    #   Reference: Paolella M 0886. Tail Estimation and Conditional Modeling
    #   of Heteroscedastic Time!Series. PhD thesis.
    #   Institute of Statistics and Econometrics. Christian Albrechts University at Kiel
    #   Parameters: mean in R; sd > 0; nu > 0; d > 0; xi > 0;

    # FUNCTION:

    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }

    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")

    # Compute auxiliary variables:
    z = (x - mean ) / sd
    n = length(z)
    arg = z
    indexLessThanZero = which (z < 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)

    # Compute the density points according to their sign
    # all b coefficients are >= 0
    if(sizeIndex == 0) {
      arg = arg / xi
      # all b coefficients are < 0
    } else if (sizeIndex == n) {
      arg = -arg * xi
      # default case. we have both pos. and neg. values
    } else if (TRUE) {
      arg[indexLessThanZero] = -arg[indexLessThanZero] * xi
      arg[-indexLessThanZero] = arg[-indexLessThanZero] / xi
    }

    # Compute density points
    k = ( ( xi + 1/xi ) * 1/d * nu^(1/d) * beta (1/d,nu) )^(-1)
    result = ( k * (1 + (arg^d) / nu )^( -nu-1/d) ) / sd
    # Log:
    if(log) result = log(result)

    # Return Value
    result
  }





# ------------------------------------------------------------------------------

pgat <-
  function(q, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {
    # A function imlemented by Thiago Sousa

    # Description:
    #   Compute the distribution for the
    #   so called t3-distribution (now called GAt)
    #   Parameters: mean in R; sd > 0; nu > 0; d > 0; xi > 0;

    # FUNCTION:

    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }

    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")

    # Define auxiliary functions
    L <- function (z, nu = nu, d = d, xi = xi)
    {
      nu / ( nu + (-z*xi)^d ) # z must be negative ( <= 0 ), but we do not check it here
    }
    U <- function (z, nu = nu, d = d, xi = xi)
    {
      pw = (z/xi)^d  # z must be negative ( <= 0 ), but we do not check it here
      return ( replace(pw / ( nu + pw ),which (pw == Inf, arr.ind = TRUE),1) )
    }

    # Compute auxiliary variables:
    z = (q - mean ) / sd
    n = length(z)
    arg = z
    indexLessThanZero = which (z <= 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)

    # Compute distribution points according to their sign
    if(sizeIndex == 0) {
      arg = 1/(1 + xi^2 ) + 1/(1 + xi^(-2) ) *
        pbeta ( U (z = arg, nu = nu, d = d, xi = xi), 1/d, nu)
    } else if (sizeIndex == n) {
      arg = 1/(1 + xi^2 ) *
        pbeta ( L (z = arg, nu = nu, d = d, xi = xi), nu, 1/d)
    } else if (TRUE) {
      arg[indexLessThanZero] = 1/(1 + xi^2 ) *
        pbeta ( L (z = arg[indexLessThanZero], nu = nu, d = d, xi = xi), nu, 1/d)
      arg[-indexLessThanZero] = 1/(1 + xi^2 ) + 1/(1 + xi^(-2) ) *
        pbeta ( U (z = arg[-indexLessThanZero], nu = nu, d = d, xi = xi), 1/d, nu)
    }

    # Return Value
    arg
  }


#------------------------------------------------------------------------------


qgat <-
  function(p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {

    # Description:
    #   Compute the quantiles for the
    #   generalized error distribution using the
    #   formula for the distribution function and
    #   the quantile function of the Beta Distribution
    #   already available in R

    # FUNCTION:

    # Define auxiliary functions
    Lp <- function (p = p, nu = nu, d = d, xi = xi)
    {
      qbeta( ( 1 + xi^2 ) * p, nu, 1/d)
    }
    Up <- function (p = p, nu = nu, d = d, xi = xi)
    {
      qbeta( ( p - 1 / ( 1 + xi^2 ) ) * ( 1 + xi^(-2) ), 1/d, nu)
    }

    # Compute quantiles located at (-Inf,0] and at (0,+Inf)
    F0 = pgat(0, mean = 0, sd = 1, nu = nu, d = d, xi = xi)
    n = length(p)
    result = rep(NA,n)
    indexLessThanF0 = which (p <= F0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanF0)
    if(sizeIndex == 0) {

      U = Up (p = p, nu = nu, d = d, xi = xi)
      result = ( U * nu / (1 - U) )^( 1/d ) * xi

    } else if (sizeIndex == n) {

      L = Lp (p = p, nu = nu, d = d, xi = xi)
      result = - ( nu/L - nu )^( 1/d ) * 1/xi

    } else if (TRUE) {

      L = Lp (p = p[indexLessThanF0], nu = nu, d = d, xi = xi)
      U = Up (p = p[-indexLessThanF0], nu = nu, d = d, xi = xi)
      result[indexLessThanF0] = - ( nu/L - nu )^( 1/d ) * 1/xi
      result[-indexLessThanF0] = ( U * nu / (1 - U) )^( 1/d ) * xi
    }

    # Return Value:
    result * sd + mean
  }







# ------------------------------------------------------------------------------


rgat <-
  function(n, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {

    # Description:
    #   Generate GAt Random values
    #   using the inverse of the distribution function.

    # FUNCTION:

    randomUnif = runif(n = n, min = 0, max = 1)
    result = qgat(p = randomUnif, mean = mean, sd = sd, nu = nu, d = d, xi = xi)

    # Return Value:
    result
  }


################################################################################


# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


##################################################################################
#  FUNCTION:             Skewed version (Fernandez, C. and Steel, M. F. J. (1998))
#                        of the standard t-Student (std) distribution
#
#  dskstd                Density for the skstd
#  pskstd                Probability function for the skstd
#  qskstd                Quantile function for the skstd
#  rskstd                Random Number Generator for the skstd
##################################################################################


dskstd <-
  function(x, mean = 0, sd = 1, nu = 3, xi = 1, log = FALSE)
  {
    # A function implemented by Thiago Sousa

    # DESCRIPTION:
    #   Compute the density for the
    #   skewed version of the sandard t-Student.
    #   The skewed t-Student distributio is already defined on
    #   the skewt R package, but it uses th non-standardized version
    #   of the t-Student distribution, i.e., the variance is df/(df-2)
    #   and it does not allow the user to specify the location or
    #   scale parameter. We are defining this slightly diffeernt
    #   version just because it is widelly used in the GARCH context.
    #   Reference: Fernandez, C. and Steel, M. F. J. (1998).
    #   On Bayesian modeling of fat tails and skewness,
    # J. Am. Statist. Assoc. 93, 359-371.

    # INPUT PARAMETERS:
    #  x - the points in which the function will be evaluated
    #  mean - the location parameter
    #  sd - the scale parameter
    #  nu - the shape parameter ( degrees of freedom) ( > 0)
    #  xi - the asymmetry parameter ( > 0 ). xi = 1 gives a symmetric distr.
    #  log - logical indicating whether or not to report the log of the density

    # FUNCTION:

    # Params:
    if (length(mean) == 4) {
      xi = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }

    # Error treatment of input parameters
    if(sd <= 0  || nu <= 2 || xi <= 0)
      stop("Failed to verify condition:
           sd <= 0  || nu <= 0 || xi <= 0")

    # Compute density points:
    z = (x - mean ) / sd
    result = sqrt( nu / ( nu - 2 ) ) *
      dskt(z * sqrt( nu / ( nu - 2 ) ), df = nu, gamma = xi)

    # Log:
    if(log) result = log(result)

    # Return Value
    result
  }



# ------------------------------------------------------------------------------


pskstd <-
  function(q, mean = 0, sd = 1, nu = 3, xi = 1)
  {
    # A function implemented by Thiago Sousa

    # DESCRIPTION:
    #   Compute the probability function for the
    #   skewed version of the sandard t-Student.

    # INPUT PARAMETERS:
    #  q - the points in which the function will be evaluated
    #  mean - the location parameter
    #  sd - the scale parameter
    #  nu - the shape parameter ( degrees of freedom) ( > 0)
    #  xi - the asymmetry parameter ( > 0 ). xi = 1 gives a symmetric distr.
    #  log - logical indicating whether or not to report the log of the density

    # FUNCTION:

    # Params:
    if (length(mean) == 4) {
      xi = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }

    # Error treatment of input parameters
    if(sd <= 0  || nu <= 2 || xi <= 0)
      stop("Failed to verify condition:
           sd <= 0  || nu <= 0 || xi <= 0")

    # Compute probability points:
    z = (q - mean ) / sd
    result = pskt( z * sqrt( nu / ( nu - 2 ) ), df = nu, gamma = xi)

    # Return Value
    result
  }


# ------------------------------------------------------------------------------


qskstd <-
  function(p, mean = 0, sd = 1, nu = 3, xi = 1)
  {
    # A function implemented by Thiago Sousa

    # DESCRIPTION:
    #   Compute the quantile function for the
    #   skewed version of the sandard t-Student.

    # INPUT PARAMETERS:
    #  p - the probabilities in which the function will be evaluated
    #  mean - the location parameter
    #  sd - the scale parameter
    #  nu - the shape parameter ( degrees of freedom) ( > 0)
    #  xi - the asymmetry parameter ( > 0 ). xi = 1 gives a symmetric distr.
    #  log - logical indicating whether or not to report the log of the density

    # FUNCTION:

    # Params:
    if (length(mean) == 4) {
      xi = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }

    # Error treatment of input parameters
    if(sd <= 0  || nu <= 2 || xi <= 0)
      stop("Failed to verify condition:
           sd <= 0  || nu <= 0 || xi <= 0")

    # Compute probability points:
    result = qskt( p, df = nu, gamma = xi) /  sqrt( nu / ( nu - 2 ) )

    # Return Value:
    result * sd + mean
  }



# ------------------------------------------------------------------------------


rskstd <-
  function(n, mean = 0, sd = 1, nu = 3, xi = 1)
  {

    # Description:
    #   Generate skstd Random values
    #   using the inverse of the distribution function.

    # FUNCTION:

    randomUnif = runif(n = n, min = 0, max = 1)
    result = qskstd(p = randomUnif, mean = mean, sd = sd, nu = nu, xi = xi)

    # Return Value:
    result
  }


################################################################################







# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .armaGarchDist          Computes density values for several
#                          conditional distributions
################################################################################




# ------------------------------------------------------------------------------
# The default is to use the dstable density from package stabledist
# If the user has the 'stable' package from JP Nolan, then the density
# is computed with the dstable.quick function.
# The configuration of the .GSgarch.dstable function was done
# inside the .onAttach function, but we get some warnings concerning
# .armaGarchDist: no visible global function definition for '.GSgarch.dstable'.
# Therefore, we decided to define the function here and to decide wether
# or not to use package 'stable' only when function .armaGarchDist is called.
# This is bad for efficiency and we must find another reasonable way to define it


.GSgarch.dstable <- function(x,alpha = 1.5, beta = 0, gamma = 1,
                             delta = 0, param = 1)
{
  return(stabledist::dstable(x, alpha, beta, gamma,
                             delta, pm = param))
}



# ------------------------------------------------------------------------------



.armaGarchDist <-
  function(z, hh, shape = 1.5, skew = 0,
           cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
           TOLG = 1e-8)
  {
    # Description:
    #   Calculates the likelihood function for a vector of points (z)
    #   according to the specified distribution.
    #   REMARKS: We choose not to avoid calling this function with out of bound parameters
    #   because in these cases we simply return a huge likelihood to guide
    #   the optimization algorithm.

    # Arguments:
    #   z - vector of points to calculate the llh.
    #   h - vector of conditional variances ????? BETTER DESCRIPTION NEEDED.
    #   cond.dist - name of the conditional distribution
    #   TOLG - general tolerance for arma-garch parameters.
    #   In the beggining it was set to 1e-5

    # FUNCTION:

    # Error treatment of input parameters
    cond.dist = match.arg(cond.dist)
    if(length(z) != length(hh))
      stop("Error: Vectors 'z' and 'hh' have different length.")
    if(sum(is.na(hh)) > 0 || min(hh) == 0)
    {
      return(1e99)
    }

    # normal conditional distribution
    if(cond.dist == "norm")
      return(-sum(log(dnorm(x = z/hh)/hh)))

    # t-student conditional distribution.
    if(cond.dist == "std")
    {
      if(!(shape > 2))
        return(Inf)
      return(-sum(log(dstd(x = z/hh, nu = shape)/hh)))
    }

    # skew t-student conditional (standardized version defined in Wurtz)
    if(cond.dist == "sstd")
    {
      if(!(shape > 2) || !(skew > 0))
        return(1e99)
      return(-sum(log(dsstd(x = z/hh, nu = shape, xi = skew)/hh)))

    }

    # skew t-student from Fernandez, C. and Steel, M. F. J. (1998)
    if(cond.dist == "skstd")
    {
      if(!(shape > 2) || !(skew > 0))
        return(1e99)
      return(-sum(log(dskstd(x = z/hh, nu = shape, xi = skew)/hh)))

    }

    # GAt distribution
    if(cond.dist == "gat")
    {
      if(!(shape[1] > 0) || !(shape[2] > 0) || !(skew > 0))
        return(1e99)
      return(-sum(log(dgat(x = z/hh, nu = shape[1], d = shape[2], xi = skew)/hh)))
    }

    # GED conditional distribution.
    if(cond.dist == "ged")
    {
      if(!(shape > 0))
        return(1e99)
      return(-sum(log(dged(x = z/hh, nu = shape)/hh)))
    }

    # GEV conditional distribution
    if(cond.dist == "gev")
    {
      # to ensure good mle properties and finiteness of the variance
      # we require shape > -0.5 ( See Jondeau et al. (XXX))
      if( (shape[1] <= -0.5) || (shape[1] >= 0.5))
        return(1e99)

      # compute density normally with shape != from zero. There was some numerical problems
      # while using function dgev from package fExtremes and thus, we decided to compute the
      # density directly, asumming that the shape parameter is different from zero. Aditionally,
      # we cannot enforce the computation with the shape parameter equal to zero, since the
      # algorithm is not able to move for the optimum in most cases.

      # when the shape parameter is small
      if( abs(shape[1]) < 1e-7 )
      {
        llh <- sum(log(hh)) + sum(z/hh) + sum(exp(-z/hh))
        return(llh)
      }

      # when the shape parameter is next to zero
      sig <- hh
      y <- 1 + shape * (z/hh)
      llh <- sum(log(sig)) + sum(y^(-1/shape)) + sum(log(y))*(1/shape + 1)
      return(llh)
    }

    # stable conditional distribution
    if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
    {
      # Compute density with package 'stable' if it is available
      if(getOption('.stableIsLoaded', default = FALSE) == TRUE)
      {
        .GSgarch.dstable <- function(x,alpha = 1.5, beta = 0, gamma = 1,
                                     delta = 0, param = 1)
        {
          return(stable::dstable.quick(x, alpha, beta, gamma,
                                       delta, param))
        }
      }




      # Return Big Values if we are out of parameter space

      if( !(shape > 1) || !(shape < 2) || !(abs(skew) < 1))
        return(1e99)

      y <- z/hh
      if(sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y)))
        return(1e99)

      # Compute density either using "stable" or "stabledist" package

      if(cond.dist == "stableS0")
        result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape,
                                           beta = skew, param = 0)/hh))

      if(cond.dist == "stableS1")
        result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape,
                                           beta = skew, param = 1)/hh))

      if(cond.dist == "stableS2")
        result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape,
                                           beta = skew, param = 2)/hh))

      # Result
      return(result)
    }
  }






################################################################################




# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA






################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterArma             Filter function for Arma, Ar or Ma process
#
################################################################################


.filterArma <- function(
    data, size.data,
    m,n,
    mu, a, b)
{
  # Arguments
  # data: vector with data
  # m,n: model order arma(m,n)
  # mu,a,b: model parameters. 'a' is the autorregressive coefficient
  # and 'b' is the moving average coeficient.
  # REMARKS: This function filters, pure 'arma' models, but
  # it can deals with other models, such as 'ar' and 'ma' models.
  # All the model parameters must be specified, regardless if they
  # exist. For example:
  # for ar(m) m >= 1 make n = 1 and b = 0.
  # for ma(n) n >= 1 make m = 1 and a = 0.

  # Return
  # the residuals 'z'

  # FUNCTION:

  # error treatment of input parameters
  if( n < 1 || n < 1 || length(a) != m || length(b) != n)
    stop("One or more of these conditions were true:
           n < 1 || n < 1 || length(a) != m || length(b) != n")

  # Initial declaration of variables
  x.init <- rep(0,m)
  z.init <- rep(0,n)

  N <- size.data
  x.ZeroMean <- data - mu
  x2 = 0
  V <- c(x.init,x.ZeroMean[1:(N-1)])
  for( i in 1:m)
  {
    x1 <- -a[i]*V
    x2 = x2 +  x1[(m-(i-1)):(m+N-i)]
  }
  x2 <- x.ZeroMean + x2

  z <- filter(x2, filter = -b,
              method = "recursive", init = z.init)

  if(length(z) != N)
    stop("Error in filtering function. length(z) != N")

  # return
  return(z)
}


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterAparch          Filtering Aparch process and its particular cases
#
################################################################################

.filterAparch <- function(
    data,init.Value = NULL,
    p,q,
    mu, omega, alpha, beta, gamma, delta)
{

  # Arguments
  # data: vector with data
  # p,q: model order. aparch(p,q)
  # mu,omega,alpha,beta,gamma,delta: model parameters
  # REMARKS: This function filters, pure 'aparch' models, but
  # it can deals with other models, such as garch, arch.
  # All the model parameters must be specified, regardless if they
  # exist. For example:
  # for garch(p,q) with p,q >= 1 make gamma = 0.
  # for arch(p) p >= 1 make q = 1, beta = 0 and gamma = 0.

  # Return
  # the series 'z' and 'h'


  # input
  # data: vector with data

  # error treatment of input parameters
  if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p
      || length(beta) != q)
    stop("One or more of these conditions were true:
          p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p
          || length(beta) != q")

  # Initial declaration of variables
  pq = max(p,q)
  z = (data-mu)
  N = length(data)
  Mean.z = mean(abs(z)^delta)

  # initializing the time series
  if(is.null(init.Value)) {
    edelta.init <- rep(Mean.z,p)
    h.init <- rep(Mean.z,q)
  } else {
    edelta.init <- rep(init.Value,p)
    h.init <- rep(init.Value,q)
  }

  edeltat = 0
  for( i in 1:p)
  {
    edelta <- alpha[i]*(c(edelta.init,((abs(z)-gamma[i]*z)^delta)[1:(N-1)]))
    edeltat = edeltat +  edelta[(p-(i-1)):(p+N-i)]
  }
  edeltat = omega + edeltat

  h <- filter(edeltat, filter = beta,
              method = "recursive", init = h.init)

  if(length(z) != length(h))
    stop("Error in filtering function. length(z) != length(h)")
  hh <- (abs(h))^(1/delta)

  # return
  cbind(z,hh)
}

################################################################################









################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterGarch11Fit      Filter function for Garch(1,1) process. This code
#                          snipet was taken from the Wurtz et al. (2006) paper
#
################################################################################

.filterGarch11Fit <- function(data, parm)
{
  mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
  z = (data-mu); Mean = mean(z^2)

  # Use Filter Representation:
  e = omega + alpha * c(Mean, z[-length(data)]^2)
  h = filter(e, beta, "r", init = Mean)
  print(h[1])
  hh = sqrt(abs(h))
  cbind(z,hh)
}

################################################################################













################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterAparchForLoop   Conditional Variance filtering -
#                          for loop - Wuertz et al. (2006)
#
################################################################################

.filterAparchForLoop <- function(
    data,h.init = 0.1,
    p,q,
    mu, omega, alpha, beta, gamma, delta)
{

  # Return
  # the series 'z' and 'hh'


  # input
  # data: vector with data

  # error treatment of input parameters
  if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p
      || length(beta) != q)
    stop("One or more of these conditions were true:
            p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p
            || length(beta) != q")

  # Initial declaration of variables
  pq = max(p,q)

  z = (data-mu)
  N = length(data)
  Mean.z = mean(abs(z)^delta)
  h = rep(h.init, pq)

  # Calculate h[(pq+1):N] recursively
  for (i in (pq+1):N )
  {
    ed = 0
    for (j in 1:p)
    {
      ed = ed+alpha[j]*(abs(z[i-j])-gamma[j]*z[i-j])^delta
    }
    h[i] = omega + ed + sum(beta*h[i-(1:q)])
  }
  if(length(z) != length(h))
    stop("Error in filtering function. length(z) != length(h)")

  hh <- (abs(h))^(1/delta)

  # return
  cbind(z,hh)
}

################################################################################



# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  gsFit                   Fits ARMA-GARCH or ARMA-APARCH model
################################################################################

gsFit <-
  function(
    formula = ~ garch(1,1),
    data,
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
    include.mean = TRUE,
    algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"),
    control = NULL,
    tolerance = NULL,
    title = NULL,
    description = NULL)
  {
    # Description:
    #     This functions reads the univariate time series and fits
    #     a ARMA-GARCH/APARCH model with conditional GEV and stable
    #     distribution.
    #     TIME SERIES MODEL:
    #         According to Wurtz et al. (2006)
    #         xt = mu + a(B)xt + b(B)et,
    #         et = zt * sigmat
    #         zt ~ Dv(0, 1)
    #         sigmat^2 = omega + sum(alphai(e(t-i))^2) + sum(betaj*sigmat-j^2)
    #     REMARKS:
    #     et for ARMA(m,n), n > 1 will be initiated as 0, i.e, e[-n+1:0] = 0.
    #     zt will be initiated as 0.
    #     ht will be initiated as 0.1 (see eq. (22) of Wurtz, 2006).
    #     VARIABLE NOTATION USED INSIDE THIS FUNCTION:
    #         N: sample size
    #         m,n: ARMA order
    #         p,q, pq = max(p,q): GARCH order. They were called u,v, uv in Wurtz et al.(2006)
    #         mu, a, b: ARMA parameters
    #         omega, alpha, gamma, beta: APARCH vector parameters
    #         garchLLH: Log Likelihood for the ARMA(m,n)-GARCH(p,q) model
    #         sigma: The scale parameter in a pure ARMA(m,n) model with innovations
    #         D(shape,skew,sigma,location = 0).
    #         llh: The negative of the log-likelihood
    #         parm: ARMA-GARCH parameters concatenated as c(mu,a,b,omega,alpha,gamma,beta)
    #         h: conditional variance. It was referenced as sigmat in model equations
    #         x: data set. It is the time series that will be modelled as ARMA-GARCH
    #         xi: GEV shape parameter. xi > -0.5 (llh) and < 0.5 (finiteness of variance, see Zhao et al. (2011))
    #         AR, MA: if AR (or MA) equals TRUE, then the model has AR (or MA) order equals to zero.
    #         Stable distributions can be specified in all parametrizations: S0, S1 and S2.
    #         Even for stable distribution llh, there's a problem for finding the estimators for
    #         alpha near to 2 and beta > 0. The bad performance on the ARMA-GARCH model was
    #         similar to the single llh estimation of i.i.d samples of stable distribution.
    #         In our simulations, the param = 2 was chosen to be the ones that performs
    #         better in these situations.
    #     IMPORTANT DETAILS FOR USING THIS FUNCTION:
    #         Some care must be taken when using this function to estimate the parameters
    #         of some models.
    #         For pure GARCH(p,q) with p,q >= 1 we need to make sure that the 'gamma' variable
    #         is a vector of length 'p' and with all entries equal to zero.
    #         For pure ARCH(p) = GARCH(p,0) we set the variable GARCH equal to TRUE to indicate
    #         that the model has order q = 0. Then, we make q = 1 to estimate the parameters
    #         of a GARCH(p,1) with beta = 0 and gamma = 0.

    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification
    #   data - vector of data
    #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
    #   include.mean - a logical, should the mean value be estimated ?
    #   algorithm - the selected algorithm
    #   cond.dist - name of the conditional distribution
    #   title - a character string which allows for a project title
    #   description - a character string which allows for a project description

    # Return:
    #   result - An object of class GEVSTABLEGARCH

    # FUNCTION:

    # Error Treatment on input parameters
    DEBUG = FALSE
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    if (!is.numeric(data) || any(is.na(data)) || any(is.null(data)) || any(!is.finite(data)))
      stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")

    # Call:
    CALL = match.call()

    # Configuring Tolerance
    if( is.null (tolerance) )
      tolerance = list( TOLG = 1e-8, TOLSTABLE = 1e-2, TOLSTATIONARITY = 1e-3 )
    if( tolerance$TOLSTATIONARITY > 0.05)
      stop("TOLSTATIONARITY can not be > 0.05")

    if(is.null(title))
      title = "ARMA-GARCH modelling"

    # Getting order model from object formula
    formula.input <- formula
    formula <- .getFormula(formula)
    m <- formula$formula.order[1]
    n <- formula$formula.order[2]
    p <- formula$formula.order[3]
    q <- formula$formula.order[4]
    APARCH <- formula$isAPARCH
    formula.mean <- ""
    formula.var <- ""
    if(m > 0 || n > 0)
      formula.mean <- formula$formula.mean
    else
      formula.mean <- ""
    if(p > 0)
      formula.var <- formula$formula.var
    else
      formula.var <- ""

    # Stop if the algorithm is not supported
    if(algorithm == "sqp.restriction")
    {
      if(formula.var == "")
        stop("sqp.restriction should only be used with GARCH/APARCH models.")
      if( any ( cond.dist == c("stableS0", "stableS2") ) )
        stop("sqp.restriction algorithm can only be used with \n stable distribution in S1 parametrization, i.e., cond.dist = 'stableS1'.")
    }

    # Stop if the user specified a pure arma model
    if(formula.var == "")
      stop("Pure ARMA model not allowed")


    # Configuring model order
    printRes = TRUE
    ARMAonly <- FALSE
    AR <- FALSE
    MA <- FALSE
    GARCH <- FALSE
    if( m == 0) AR <- TRUE
    if( n == 0) MA <- TRUE
    if( q == 0) GARCH <- TRUE
    optim.finished <- FALSE
    if (AR == TRUE)
      m <- 1
    if( MA == TRUE)
      n <- 1
    if( (p == 0) && (q == 0))
      ARMAonly = TRUE
    if (GARCH == TRUE && !ARMAonly)
      q <- 1

    # Initial configurations
    data <- data;
    N <- length(data)
    out <- NULL # output of the gsFit function
    messages <- NULL # important messages collected during estimation
    out$order <- c(m,n,p,q,include.mean,APARCH)
    TMPvector <- c(if(m != 0 || n != 0) c("arma(",m,",",n,")-"),
                   if(APARCH==FALSE)c("garch(",p,",",q,")"),
                   if(APARCH==TRUE)c("aparch(",p,",",q,")"))
    TMPorder <- paste(TMPvector, collapse="")
    TMPvectorintercept <- c("include.mean:",if(include.mean==TRUE)"TRUE",
                            if(include.mean==FALSE)"FALSE")
    TMPintercept <- paste(TMPvectorintercept, collapse="")
    out$model <- paste(TMPorder,"##",TMPintercept, collapse="")
    out$cond.dist <- cond.dist
    out$data <- data
    optim.finished <- FALSE
    if(cond.dist == "gat")
      lengthShape = 2
    else
      lengthShape = 1

    if(DEBUG)
      print(c("lengthShape",lengthShape))
    # This function checks if the model is an stationary ARMA process.
    arCheck <- function(ar) {
      p <- max(which(c(1, -ar) != 0)) - 1
      if (!p)
        return(TRUE)
      all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
    }


    # BEGIN: garch Likelihood ARMA-APARCH or pure GARCH process
    ###########################################################
    garchLLH = function(parm){

      # check if some parameters are NAN
      if(sum(is.nan(parm)) != 0) {return(1e99)}

      # Getting parameters from parm vector
      mu <- parm[1];
      a <- parm[(1+1):(2+m-1)]
      b <- parm[(1+m+1):(2+m+n-1)]
      omega <- parm[1+m+n+1]
      alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
      gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
      beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
      delta <- parm[2+m+n+2*p+q+1]
      skew <- parm[3+m+n+2*p+q+1]
      shape <- parm[(4+m+n+2*p+q+1):(4+m+n+2*p+q+lengthShape)]


      # Configure 'shape' for the 'GAt' distribution
      #if(cond.dist == "gat")
      #    shape <- parm[4+m+n+2*p+q+1,4+m+n+2*p+q+1+1]

      # Configuring delta and gamma for Garch estimation
      if( !APARCH)
      {
        gm = rep(0,p);
        delta = 2;
        if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
          delta = 1
      }

      # Setting parameters to accept tapper off MA, AR or GARCH coefficients
      if( AR == TRUE)
        a <- 0
      if( MA == TRUE)
        b <- 0
      if( GARCH == TRUE)
        beta <- 0
      if (include.mean == FALSE)
        mu <- 0

      # Avoid being out of parameter space
      cond.general <- FALSE
      cond.normal <- FALSE
      cond.student <- FALSE
      cond.gev <- FALSE
      cond.stable <- FALSE
      parset <- c(omega,alpha,if(!GARCH) beta,delta)
      cond.general <- any(parset < tolerance$TOLG)

      if (cond.general || cond.student || cond.gev || cond.stable)
      {
        return(1e99)
      }

      # ARMA stationarity condition check
      if(!arCheck(a))
        return(1e99)


      # Filters the Time series to obtain the i.i.d. sequence of
      # 'innovations' to evaluate the Log-likelihood function
      if(AR == TRUE && MA  == TRUE)
      {
        filteredSeries <- .filterAparch(data = data,p = p,q = q,
                                        mu = mu, omega = omega, alpha = alpha, beta = beta, gamma = gm, delta = delta)
        z <- filteredSeries[,1]
        hh <- filteredSeries[,2]
      }
      if(AR == FALSE || MA == FALSE)
      {
        filteredArma <- .filterArma(data = data, size.data = N, m = m, n = n, mu = mu, a = a, b = b)
        filteredSeries <- .filterAparch(data = filteredArma,p = p,q = q,
                                        mu = 0, omega = omega, alpha = alpha, beta = beta, gamma = gm, delta = delta)
        z <- filteredSeries[,1]
        hh <- filteredSeries[,2]

      }


      # get output Residuals
      if (optim.finished)
      {
        out$residuals <<- as.numeric(z)
        out$sigma.t <<- as.numeric(hh)
        out$h.t <<- as.numeric(hh^delta)
      }

      # Return llh function
      llh.dens <- .armaGarchDist(z = z, hh = hh, shape = shape,
                                 skew = skew, cond.dist = cond.dist)
      llh <- llh.dens
      if (is.nan(llh) || is.infinite(llh) || is.na(llh))
      {
        llh <- 1e99
      }
      llh
    }
    # END: garch Likelihood ARMA-APARCH or pure GARCH process
    ###########################################################


    # Getting start, lower and upper bound for parameters to perform optimization
    start <- .getStart(data = data,m = m,n = n,p = p,q = q,AR = AR,
                       MA = MA, cond.dist = cond.dist,
                       TOLG = tolerance$TOLG, TOLSTABLE = tolerance$TOLSTABLE)
    if(DEBUG)
    {
      print("start")
      print(start)
    }

    # Function that evaluate stationarity conditions to guide parameter estimation.
    garch.stationarity <- function(parm)
    {

      # Getting parameters from parm vector
      # mu <- parm[1];
      # a <- parm[(1+1):(2+m-1)]
      # b <- parm[(1+m+1):(2+m+n-1)]
      omega <- parm[1+m+n+1]
      alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
      gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
      beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
      delta <- parm[2+m+n+2*p+q+1]
      skew <- parm[3+m+n+2*p+q+1]
      shape <- parm[(4+m+n+2*p+q+1):(4+m+n+2*p+q+lengthShape)]

      model = list(omega = omega, alpha = alpha, gm = gm, beta = beta,
                   delta = delta, skew = skew, shape = shape)


      # Return
      .stationarityAparch(model = model,
                          formula = formula,
                          cond.dist = cond.dist)
    }

    # Optimization procedure using selected algorithms
    modelLLH <- garchLLH


    if (algorithm == "nlminb")
    {
      fit1 <- nlminb(start[1,], objective = modelLLH,
                     lower = start[2,], upper = start[3,],
                     control = control)
      out$llh <- fit1$objective
      out$par <- fit1$par
      out$hessian <- fit1$hessian
      out$convergence <- fit1$convergence
    }
    if (algorithm == "nlminb+nm")
    {
      fit1.partial <- nlminb(start[1,], objective = modelLLH,
                             lower = start[2,], upper = start[3,],
                             control = control)
      fit1 <- optim(par = fit1.partial$par, fn = modelLLH,
                    method = "Nelder-Mead", hessian = TRUE)
      out$llh <- fit1$value
      out$par <- fit1$par
      out$hessian <- fit1$hessian
      out$convergence <- fit1$convergence
    }
    if (algorithm == "sqp")
    {
      fit1 <- solnp(pars = start[1,], fun = modelLLH,
                    LB = start[2,], UB = start[3,], control = control)
      out$llh <- fit1$values[length(fit1$values)]
      out$par <- fit1$pars
      out$hessian <- fit1$hessian
      out$convergence <- fit1$convergence
    }
    if (algorithm == "sqp.restriction")
    {
      fit1 <- solnp(pars = start[1,], fun = modelLLH, ineqfun = garch.stationarity, ineqLB = 0,
                    ineqUB = 1-tolerance$TOLSTATIONARITY, LB = start[2,], UB = start[3,],
                    control = control)
      out$llh <- fit1$values[length(fit1$values)]
      out$par <- fit1$pars
      out$hessian <- fit1$hessian
      sizeHessian = length(out$hessian[1,])-1
      out$hessian = out$hessian[1:sizeHessian,1:sizeHessian]
      out$convergence <- fit1$convergence
    }
    if ( (any(c("sqp","nlminb")  == algorithm)) && !is.numeric(try(sqrt(diag(solve(out$hessian))), silent = TRUE)))
    {
      # Call the Nelder-Mead method just to calculate the hessian matrix
      # This is due to the fact that the hessian matrix returned
      # by the sqp routine is often difficult to invert, which does not
      # happen with the hessian returned by the optim routine
      # using the "Nelder-Mead" algorithm. Note that
      # we call it with only one interation, so that the algorithm does
      # not change the parameter values and only calculate the hessian.
      fit1.partial <- optim(par = fit1$par, fn = modelLLH,
                            method = "Nelder-Mead", hessian = TRUE,
                            control = list(maxit = 1))

      # Check if the Nelder-Mead method changed our estimated parameters
      if(sum(abs(fit1.partial$par-fit1$par)) == 0)
      {
        out$hessian = fit1.partial$hessian
        if(DEBUG)
          print(abs(fit1.partial$par-fit1$par))
      }
    }
    if ( (algorithm == "sqp.restriction") )
    {
      # Call the Nelder-Mead method just to calculate the hessian matrix
      # This is due to the fact that the hessian matrix returned
      # by the sqp routine is often difficult to invert, which does not
      # happen with the hessian returned by the optim routine
      # using the "Nelder-Mead" algorithm. Note that
      # we call it with only one interation, so that the algorithm does
      # not change the parameter values and only calculate the hessian.
      fit1.partial <- optim(par = fit1$par, fn = modelLLH,
                            method = "Nelder-Mead", hessian = TRUE,
                            control = list(maxit = 1))

      # Check if the Nelder-Mead method changed our estimated parameters
      if(sum(abs(fit1.partial$par-fit1$par)) == 0)
      {
        out$hessian = fit1.partial$hessian
        if(DEBUG)
          print(abs(fit1.partial$par-fit1$par))
      }
    }

    if(DEBUG)
      print(fit1)

    # Configuring the convergence variable
    if(out$llh == 1e99)
      out$convergence = 1

    # Organizing the output of the program
    optim.finished = TRUE

    # Call garchLLH function to update the values of the ARMA residuals
    # and the GARCH/APARCH volatility.
    modelLLH(out$par)

    # Creating index to create a vector with the estimated parameters.
    if(!ARMAonly)
    {
      outindex <- c(if(include.mean) 1,
                    if(AR == FALSE) (1+1):(2+m-1),
                    if(MA == FALSE) (1+m+1):(2+m+n-1),
                    if(!ARMAonly) (1+m+n+1),
                    if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
                    if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
                    if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
                    if(APARCH) (2+m+n+2*p+q+1),
                    if(any(c("sstd","skstd","stableS0","stableS1","stableS2","gat")  == cond.dist)) (3+m+n+2*p+q+1),
                    if(any(c("std","gev","stableS0","stableS1","stableS2","sstd","skstd","ged","gat")  == cond.dist))
                      (4+m+n+2*p+q+1):(4+m+n+2*p+q+lengthShape))
    } else {
      outindex <- c(if(include.mean) 1,
                    if(AR == FALSE) (1+1):(2+m-1),
                    if(MA == FALSE) (1+m+1):(2+m+n-1),
                    if(!ARMAonly) (1+m+n+1),
                    if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
                    if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
                    if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
                    if(APARCH) (2+m+n+2*p+q+1),
                    if(any(c("sstd","skstd","stableS0", "stableS1", "stableS2","gat")  == cond.dist)) (1+m+n+2*p+q+1),
                    if(any(c("std","gev","stableS0", "stableS1", "stableS2","sstd","skstd","ged","gat")  == cond.dist))
                      (2+m+n+2*p+q+1):(2+m+n+2*p+q+lengthShape),
                    length(out$par))
    }

    outnames <- c(if(include.mean) "mu",
                  if( AR == FALSE) paste("ar", 1:m, sep = ""),
                  if(MA == FALSE) paste("ma", 1:n, sep = ""),
                  if(!ARMAonly) "omega",
                  if(!ARMAonly) paste("alpha", 1:p, sep = ""),
                  if(APARCH) paste("gamma", 1:p, sep = ""),
                  if(!GARCH) paste("beta", 1:q, sep = ""),
                  if(APARCH) "delta",
                  if(any(c("sstd","stableS0", "stableS1", "stableS2","gat","skstd")  == cond.dist)) "skew",
                  if(any(c("std","gev","stableS0", "stableS1", "stableS2","sstd","ged","gat","skstd")  == cond.dist))
                    paste("shape", 1:lengthShape, sep = ""),
                  if(ARMAonly) "sigma")

    if(DEBUG)
    {
      print(c("out",out))
      print(c("outindex",outindex))
    }

    out$par <- out$par[outindex]
    names(out$par) <- outnames
    out$hessian <- out$hessian[outindex,outindex]
    nParam <- length(out$par)
    out$aic  = 2*out$llh + 2*nParam # out$llh is the negative of the log-likelihood
    out$aicc = 2*out$llh + 2*(nParam+1)*N/(N - nParam - 2)
    out$bic =  2*out$llh + nParam*log(N)
    out$ics = c(out$aic,out$bic,out$aicc)
    names(out$ics) <- c("AIC","BIC","AICc")

    # Print Summary
    if ( printRes )
    {
      solveHessianFailed = FALSE
      out$se.coef <- 0
      out$se.coef <- try(sqrt(diag(solve(out$hessian))), silent = TRUE)
      if(!is.numeric(out$se.coef))
      {
        solveHessianFailed = TRUE
        messages$solving.hessian.matrix = "Error inverting Hessian Matrix"
        out$matcoef <- cbind(out$par, rep(NA,length(out$par)), rep(NA,length(out$par)), rep(NA,length(out$par)))
        dimnames(out$matcoef) = dimnames(out$matcoef) =
          list(names(out$par), c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
      }
      else
      {
        out$tval <- try(out$par/out$se.coef, silent = TRUE)
        out$matcoef = cbind(out$par, if(is.numeric(out$se.coef)) out$se.coef, if(is.numeric(out$tval)) out$tval,
                            if(is.numeric(out$tval)) 2*(1-pnorm(abs(out$tval))))
        dimnames(out$matcoef) = list(names(out$tval),
                                     c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
      }
      cat("\nFinal Estimate of the Negative LLH:\n")
      cat("-LLH:",out$llh)

      if(out$convergence == 0)
        messages$optimization.algorithm = "Algorithm achieved convergence"
      else
        messages$optimization.algorithm = "Algorithm did not achieved convergence"

      cat("\nCoefficient(s):\n")
      printCoefmat(round(out$matcoef,digits=6), digits = 6, signif.stars = TRUE)
    }

    out$order <- c(formula$formula.order[1],formula$formula.order[2],
                   formula$formula.order[3],formula$formula.order[4])
    names(out$order) <- c("m","n","p","q")
    fit <- list(par = out$par, llh = out$llh, hessian = out$hessian, ics = out$ics,
                order = out$order, cond.dist = cond.dist, se.coef = out$se.coef,
                tval = out$tval, matcoef = out$matcoef)

    # creating the output object as in fGarch package...
    new("GEVSTABLEGARCH", call = as.call(match.call()), formula = formula.input,
        method = "Max Log-Likelihood Estimation",
        convergence = out$convergence,
        messages = messages,
        data = data, fit = fit, residuals = out$residuals,
        h.t = out$h.t, sigma.t = as.vector(out$sigma.t), title = as.character(title),
        description = as.character(description))
  }

# -----------------------------------------------------------------------------




GSgarch.dstable <<- function(x,alpha = 1.5, beta = 0, gamma = 1,
                             delta = 0, param = 1)
{
  return(stabledist::dstable(x, alpha, beta, gamma,
                             delta, pm = param))
}
if(getOption('.stableIsLoaded', default = FALSE) == TRUE)
{
  GSgarch.dstable <<- function(x,alpha = 1.5, beta = 0, gamma = 1,
                               delta = 0, param = 1)
  {
    return(stable::dstable.quick(x, alpha, beta, gamma,
                                 delta, param))
  }
}



# -----------------------------------------------------------------------------




# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA





################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .getFormula             Gets the formula.mean and formula.variance from
#                          the object formula.
################################################################################


.getFormula <-
  function(
    formula)
  {
    # Description:
    #   This functions reads formula object and converts it into a list
    #   containing the separeted mean and variance arguments.
    #   Examples from output:
    #     ~arma(1,1)+garch(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(1,1)
    #     ~aparch(1,1): formula.mean = ~arma(0,0); formula.variance = ~aparch(1,1)
    #     ~arch(1): formula.mean = ~arma(0,0); formula.variance = ~arch(1)
    #     ~arma(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
    #     ~ar(1): formula.mean = ~ar(1); formula.variance = ~garch(0,0)
    #     ~ma(1): formula.mean = ~ma(1); formula.variance = ~garch(0,0)

    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification

    # Return:
    #   A list containing two elements, formula.mean and formula.variance,
    #   formula.order and the boolean formula.isAPARCH

    # FUNCTION:

    # Initial variable declaration
    allLabels = attr(terms(formula), "term.labels")
    formulas.mean.allowed = c("arma")
    formulas.variance.allowed = c("garch","aparch")
    formulaOK <- TRUE
    checkFormulaMean <- ""
    checkFormulaVariance <- ""
    isAPARCH = FALSE

    # Error treatment of input parameters
    if( (length(allLabels) != 1 ) && (length(allLabels) != 2 ) )
      formulaOK <- FALSE

    # Formula of type: ~ formula1 + formula2
    else if (length(allLabels) == 2)
    {
      formula.mean = as.formula(paste("~", allLabels[1]))
      formula.var = as.formula(paste("~", allLabels[2]))
      checkFormulaMean = rev(all.names(formula.mean))[1]
      checkFormulaVariance = rev(all.names(formula.var))[1]
      if( !any(formulas.mean.allowed == checkFormulaMean) ||
          !any(formulas.variance.allowed == checkFormulaVariance))
        formulaOK <- FALSE
    }
    # Formula of type: ~formula1
    else if (length(allLabels) == 1)
    {

      # pure 'garch' or 'aparch'
      if(grepl("arch", attr(terms(formula), "term.labels")))
      {
        formula.mean = as.formula("~ arma(0, 0)")
        formula.var = as.formula(paste("~", allLabels[1]))
        checkFormulaVariance = rev(all.names(formula.var))[1]
        if(!any(formulas.variance.allowed == checkFormulaVariance))
          formulaOK <- FALSE
      }
      else # pure 'ar', 'ma' or 'arma' model.
      {
        formula.mean = as.formula(paste("~", allLabels[1]))
        formula.var = as.formula("~ garch(0, 0)")
        checkFormulaMean = rev(all.names(formula.mean))[1]
        if(!any(formulas.mean.allowed == checkFormulaMean))
          formulaOK <- FALSE
      }
    }

    # Check if we are fitting "aparch" model
    if(checkFormulaVariance == "aparch")
      isAPARCH = TRUE

    # Get model order and check if they were specified correctly
    if(formulaOK == TRUE)
    {
      model.order.mean =
        as.numeric(strsplit(strsplit(strsplit(as.character(formula.mean),
                                              "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
      model.order.var =
        as.numeric(strsplit(strsplit(strsplit(as.character(formula.var),
                                              "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
      if( (length(model.order.mean) != 2) || (length(model.order.mean) != 2))
        formulaOK <- FALSE
    }

    # Check if model order was specified correctly.
    if(formulaOK == TRUE)
    {
      m = model.order.mean[1]
      n = model.order.mean[2]
      p = model.order.var[1]
      q = model.order.var[2]
      if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 ||
         any (c(m,n,p,q) < 0) || (p == 0 && q != 0))
        formulaOK <- FALSE
    }


    # Stop if formula was not specified correctly
    if(formulaOK == FALSE)
      stop ("Invalid Formula especification.
            Formula mean must be 'arma' and
            Formula Variance must be one of: garch or aparch
            For example:
            ARMA(1,1)-GARCH(1,1):  ~arma(1,1)+garch(1,1),
            AR(1)-GARCH(1,1):      ~arma(1,0)+garch(1,1),
            MA(1)-APARCH(1,0):     ~arma(0,1)+aparch(1,0),
            ARMA(1,1):             ~arma(1,1),
            ARCH(2):               ~garch(1,0),
            For more details just type: ?gsFit")

    # Return
    list(formula.mean = formula.mean,formula.var = formula.var,
         formula.order = c(m,n,p,q), isAPARCH = isAPARCH)
  }

################################################################################




# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .getOrder        Return a matrix with parameter order to be used
#                          inside function GSgarch.FitAIC
################################################################################


.getOrder <-
  function(order.max = c(1,1,1,1))
  {

    # Description:
    #   Iterates over the parameters to create a vector with parameter
    #   orders (like (1,1,0,1)) to use inside function GSgarch.FitAIC

    # Arguments:
    #   nMAX, mMAX, pMAX, qMAX - maximum order to be estimated

    # Return:
    #   arma.garch.order - A matrix with several model orders in the
    #   format [m,n,p,q]

    # FUNCTION:

    # error treatment on input parameters
    m = order.max[1]; n = order.max[2]; p = order.max[3]; q = order.max[4]
    if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 ||
       any (c(m,n,p,q) < 0) || (p == 0 && q != 0) ||
       any (c(m,n,p,q) > 10) )
      stop ("Invalid ARMA-GARCH order. We allow pure GARCH or APARCH. AR/MA/ARMA-GARCH/APARCH models.
              The order of the parameters could be set up to 10.")
    arma.garch.order <- c()
    for(i1 in 0:m)
    {
      for(i2 in 0:n)
      {
        for(i3 in 1:p)
        {
          for(i4 in 0:q)
          {
            ord <- c(i1,i2,i3,i4)
            arma.garch.order <- rbind(arma.garch.order,c(i1,i2,i3,i4))
          }
        }
      }
    }
    return(arma.garch.order)
  }


################################################################################



#Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .getStart               Returns initial and boundary values to
#                          perform optimization
################################################################################


.getStart <- function(data,m,n,p,q, AR = FALSE, MA = FALSE,
                      cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
                      TOLG = 1e-7, TOLSTABLE = 2e-2)
{

  # Description:
  #   Get initial values to start the estimation and the bounds for
  #   parameters to be estimated inside GSgarch.Fit function.
  #   Remarks: This function tunes initial parameters to perform optimization
  #   The ARMA coefficients are the ones returned by the "arima" function
  #   adjusted parameters functions from package ("arima" belongs to package "stats" from R)
  #   For GARCH(p,q) Wurtz et al. (2006) suggested
  #   using omega = 0.1, gamma = 0, alpha = 0.1/p and beta = 0.8/q
  #   delta is chosen to be initially equals to 1.5 in almost all the cases
  #   Keep in mind that delta < alpha for the stable case.
  #   The arma order passed to this function will be 1 even though the
  #   process does not possess the AR or MA component. Therefore,
  #   it is not a problem to have AR our MA input parameters equal to FALSE
  #   even thought the model order says on the contrary.

  # Arguments:
  #   data - vector of data
  #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
  #   AR - boolean value that indicates whether we have a model
  #   with th Autoregressive part included
  #   MA - boolean value that indicates whether we have a model
  #   with the Moving Average part included
  #   ARMAonly - Indicates whether we have a pure ARMA model
  #   cond.dist - name of the conditional distribution, one of
  #       gev, stable, norm, std, sstd, ged
  #   TOLSTABLE - boundary tolerance. Should be greater than GSstable.tol
  #   TOLG - pper and lower bounds tolerance. Should be greater than tol

  # Return:
  #   Result - A tree columns matrix with columns representing the start, lower and upper
  #     bounds of the parameter set

  # FUNCTION:

  # Error control of input parameters
  if (m < 0 || n < 0 || m %% 1 != 0 || n %% 1 != 0)
    stop("'m' and 'n' need to be integers greater than zero")

  if (m == 0 || n == 0)
    stop("Expects 'm' and 'n' different from zero")

  if ( (m != 1 && AR == TRUE)  || (n != 1 && MA == TRUE) )
    stop("If AR = TRUE, 'm' should be 1 and if MA = TRUE 'n' should be 1")

  if (p < 0 || q < 0 || p %% 1 != 0 || q %% 1 != 0)
    stop("'p' and 'q' need to be integers greater than zero")

  if(p == 0 && q != 0)
    stop("Invalid Garch(p,q) order")

  if( !is.numeric(data) || !is.vector(data))
    stop("data set must be a numerical one dimensional vector")

  # Initial variable declaration
  cond.dist = match.arg(cond.dist)
  cond.dist.list = c("stable", "gev", "gat", "norm", "std", "sstd", "skstd", "ged")
  Mean <- mean(data)
  Var <- var(data)
  Dispersion <- mean(abs(data-Mean))
  arima.fit <- c()
  arima.fit.try <- "empty"
  arima.m <- m
  arima.n <- n
  if(AR == TRUE) # we don't have the AR part
    arima.m <- 0
  if(MA == TRUE) # we don't have the MA part
    arima.n <- 0

  # Try arima fit function to get initial arma parameters
  try(arima.fit.try <- as.vector(arima(data,order = c(arima.m, 0, arima.n))$coef), silent = TRUE)
  if( is.numeric(arima.fit.try) )
  {
    arima.fit <- arima.fit.try
    if(AR == FALSE && MA == FALSE)
    {
      ar.init <- arima.fit[1:arima.m]
      ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
    }
    if(AR == TRUE && MA == FALSE)
    {
      ar.init <- 0
      ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
    }
    if(AR == FALSE && MA == TRUE)
    {
      ar.init <- arima.fit[1:arima.m]
      ma.init <- 0
    }
    if(AR == TRUE && MA == TRUE)
    {
      ar.init <- 0
      ma.init <- 0
    }
    mean.init <- arima.fit[arima.m+arima.n+1]
  } else {
    mean.init <- Mean
    ar.init <- rep(0,m)
    ma.init <- rep(0,n)
  }

  # START VALUES

  arma.start = c(ar.init, ma.init)
  mu.start = mean.init
  gm.start = rep(0, p)

  omega.start = list(
    "stableS0" = 0.1 * Dispersion,
    "stableS1" = 0.1 * Dispersion,
    "stableS2" = 0.1 * Dispersion,
    "gev" = 0.1 * Var,
    "gat" = 0.1 * Var,
    "norm" = 0.1 * Var,
    "std" = 0.1 * Var,
    "sstd" = 0.1 * Var,
    "skstd" = 0.1 * Var,
    "ged" = 0.1 * Var)

  alpha.start = list(
    "stableS0" = rep(0.1/p, p),
    "stableS1" = rep(0.1/p, p),
    "stableS2" = rep(0.1/p, p),
    "gev" = rep(0.05/p, p),
    "gat" = rep(0.1/p, p),
    "norm" = rep(0.1/p, p),
    "std" = rep(0.1/p, p),
    "sstd" = rep(0.1/p, p),
    "skstd" = rep(0.1/p, p),
    "ged" = rep(0.1/p, p))

  beta.start = list(
    "stableS0" = rep(0.8/q, q),
    "stableS1" = rep(0.8/q, q),
    "stableS2" = rep(0.8/q, q),
    "gev" = rep(0.7/q, q),
    "gat" = rep(0.8/q, q),
    "norm" = rep(0.8/q, q),
    "std" = rep(0.8/q, q),
    "sstd" = rep(0.8/q, q),
    "skstd" = rep(0.8/q, q),
    "ged" = rep(0.8/q, q))

  delta.start = list(
    "stableS0" = 1.05,
    "stableS1" = 1.05,
    "stableS2" = 1.05,
    "gev" = 2,
    "gat" = 2,
    "norm" = 2,
    "std" = 2,
    "sstd" = 2,
    "skstd" = 2,
    "ged" = 2)

  skew.start = list(
    "stableS0" = 0,
    "stableS1" = 0,
    "stableS2" = 0,
    "gev" = 1,
    "gat" = 1,
    "norm" = 1,
    "std" = 1,
    "sstd" = 1,
    "skstd" = 1,
    "ged" = 1)

  shape.start = list(
    "stableS0" = 1.9,
    "stableS1" = 1.9,
    "stableS2" = 1.9,
    "gev" = 0, # numerical tests showed that 0.01 is a good starting parameter.
    "gat" = c(2, 4),
    "norm" = 1,
    "std" = 4,
    "sstd" = 4,
    "skstd" = 4,
    "ged" = 4)

  # LOWER BOUNDS

  mu.lower = min(- 100 * abs ( mu.start ), -100)
  arma.lower = rep ( - 100, m + n )
  omega.lower = TOLG
  alpha.lower = rep ( TOLG, p )
  beta.lower = rep ( TOLG, q )
  gm.lower = rep( - 1 + TOLG, p)

  delta.lower = list(
    "stableS0" = 1,
    "stableS1" = 1,
    "stableS2" = 1,
    "gev" = TOLG,
    "gat" = TOLG,
    "norm" = TOLG,
    "std" = TOLG,
    "sstd" = TOLG,
    "skstd" = TOLG,
    "ged" = TOLG)

  skew.lower = list(
    "stableS0" = - 1 + TOLSTABLE,
    "stableS1" = - 1 + TOLSTABLE,
    "stableS2" = - 1 + TOLSTABLE,
    "gev" = 0,
    "gat" = TOLG,
    "norm" = 0,
    "std" = 0,
    "sstd" = TOLG,
    "skstd" = TOLG,
    "ged" = 0)

  shape.lower = list(
    "stableS0" = 1 + TOLSTABLE,
    "stableS1" = 1 + TOLSTABLE,
    "stableS2" = 1 + TOLSTABLE,
    "gev" = - 0.5 + TOLG, # to ensure good MLE properties. See Jondeau et al.
    "gat" = c ( TOLG, TOLG),
    "norm" = 0,
    "std" = 2 + TOLG,
    "sstd" = 2 + TOLG,
    "skstd" = 2 + TOLG, # to ensure finiteness of variance
    "ged" = TOLG)

  # UPPER BOUNDS

  mu.upper = max( 100 * abs ( mu.start ), 100)
  arma.upper = rep ( 100, m + n )
  omega.upper = 100 * abs ( Dispersion )
  alpha.upper = rep ( 1 - TOLG, p )
  beta.upper = rep ( 1 - TOLG, q )
  gm.upper = rep (1 - TOLG, p)

  delta.upper = list(
    "stableS0" = 2 - TOLSTABLE,
    "stableS1" = 2 - TOLSTABLE,
    "stableS2" = 2 - TOLSTABLE,
    "gev" = 100,
    "gat" = 100,
    "norm" = 100,
    "std" = 100,
    "sstd" = 100,
    "skstd" = 100,
    "ged" = 100)

  skew.upper = list(
    "stableS0" = 1 - TOLSTABLE,
    "stableS1" = 1 - TOLSTABLE,
    "stableS2" = 1 - TOLSTABLE,
    "gev" = 2,
    "gat" = 100,
    "norm" = 2,
    "std" = 2,
    "sstd" = 100,
    "skstd" = 100,
    "ged" = 2)

  shape.upper = list(
    "stableS0" = 2 - TOLSTABLE,
    "stableS1" = 2 - TOLSTABLE,
    "stableS2" = 2 - TOLSTABLE,
    "gev" = 0.5 - TOLG, # to ensure finiteness of the variance and mean.
    "gat" = c ( 100, 100),
    "norm" = 2,
    "std" = 100,
    "sstd" = 100,
    "skstd" = 100,
    "ged" = 100)

  # CHECK IF THE STARTING MODELS IS STATIONARY ( ONLY FOR THE sqp.restriction ALGORITHM )
  if( ! any ( cond.dist == c("stableS0", "stableS2") ) ){

    start.model.persistency = 1
    start.model.persistency = sum( alpha.start[[cond.dist]] * gsMomentAparch(
      cond.dist = cond.dist, shape = shape.start[[cond.dist]],
      skew = skew.start[[cond.dist]], delta = delta.start[[cond.dist]], gm = 0) ) +
      sum ( beta.start[[cond.dist]] )
    # The real condition is 1 but we let the user to set the limit between 0.95 and 0.999...
    if(start.model.persistency >= 0.95)
    {
      print(start.model.persistency)
      print(paste("The starting model with conditional",cond.dist, " is not stationary."))
      stop("Change the starting value of the parameters for the reported model.")
    }
  }

  # CONSTRUCT THE RESULT

  start = c ( mu.start, arma.start, omega.start[[cond.dist]], alpha.start[[cond.dist]],
              gm.start, beta.start[[cond.dist]], delta.start[[cond.dist]],
              skew.start[[cond.dist]], shape.start[[cond.dist]])

  lower = c ( mu.lower, arma.lower, omega.lower, alpha.lower,
              gm.lower, beta.lower, delta.lower[[cond.dist]],
              skew.lower[[cond.dist]], shape.lower[[cond.dist]])

  upper = c ( mu.upper, arma.upper, omega.upper, alpha.upper,
              gm.upper, beta.upper, delta.upper[[cond.dist]],
              skew.upper[[cond.dist]], shape.upper[[cond.dist]])

  namesStart = c("mu", paste("ar", 1:m, sep = ""),
                 paste("ma", 1:n, sep = ""),
                 "omega", paste("alpha", 1:p, sep = ""),
                 paste("gm", 1:p, sep = ""),
                 paste("beta", 1:q, sep = ""), "delta","skew",
                 paste("shape", 1:length(shape.start[[cond.dist]]), sep = ""))

  # Create result
  result = rbind(start,lower,upper)
  colnames(result) = namesStart

  # Return
  result
}



# ------------------------------------------------------------------------------


################################################################################



# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .onAttach               Initial configurations for stable distribution
#                          computation
################################################################################


########################################################
# Initial configurations for stable distribution
# computation
########################################################


.onAttach <-
  function(libname, pkgname)
  {
    options(.stableIsLoaded=FALSE)
    if(length(find.package("stable",quiet = TRUE)))
      options(.stableIsLoaded=TRUE)
  }


################################################################################


# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  gsSelect                Find best fitted model according to AIC criterion
################################################################################

gsSelect <-
  function(
    data,
    order.max = c(1,1,1,1),
    selection.criteria = c("AIC", "AICc", "BIC"),
    is.aparch = FALSE,
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
    include.mean = TRUE,
    algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"),
    ...)
  {


    # TESTING
    #data = x; order.max = c(1,1,1,1); is.aparch = TRUE; cond.dist = "norm"; algorithm = "sqp"
    # Description:
    #   Estimate several models using GSgarch.Fit function and evaluates
    #   the AIC to decide which one is the best model

    # Arguments:
    #   data - vector of data
    #   intercept - a logical, should the mean value be estimated ?
    #   nMAX, mMAX, pMAX, qMAX - maximum order to be estimated
    #   for the model ARMA(m,n)-GARCH/APARCH(p,q)
    #   cond.dist - name of the conditional distribution, one of
    #       gev, stable, norm, std, sstd
    #   algorithm -
    #   APARCH - boolean, indicates whether the model is an APARCH model or not

    # Return:
    #   fit.min - The best model return from function GSgarch.Fit

    # FUNCTION:

    # error treatment on input parameters
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    selection.criteria = match.arg(selection.criteria)
    if( !is.numeric(data) || !is.vector(data))
      stop("data set must be a numerical one dimensional vector")

    # begin of function
    goodness.of.fit.min <- 1e99
    goodness.of.fit.current <- 1e99
    fit.min <- list()
    fit <- list()
    order.list <- .getOrder(order.max = order.max)
    order.list.size <- length(order.list[,1])
    for( i in 1:order.list.size)
    {
      m = order.list[i,1]; n = order.list[i,2]; p = order.list[i,3]; q = order.list[i,4]

      if(is.aparch == TRUE)
      {
        formula = as.formula(paste ("~ arma(",m,", ",n,") + aparch(",p,", ",q,")", sep = "",
                                    collapse = NULL))
      }
      else {
        formula = as.formula(paste ("~ arma(",m,", ",n,") + garch(",p,", ",q,")", sep = "",
                                    collapse = NULL))
      }
      cat("\n------------------------------------------------------------------------------------------\n")
      cat(paste("Model: ", paste(formula)[2], " with '",
                cond.dist, "' conditional distribution",sep = ""))
      cat("\n------------------------------------------------------------------------------------------\n")
      fit = gsFit (data = data, formula = formula, cond.dist = cond.dist,
                   algorithm = algorithm, include.mean = include.mean, ...)

      goodness.of.fit.current = fit@fit$ics[selection.criteria]
      if (goodness.of.fit.current < goodness.of.fit.min)
      {
        fit.min <- fit
        goodness.of.fit.min <- fit@fit$ics[selection.criteria]
      }
    }
    cat("\n------------------------------------------------------------------------------------------\n")
    cat(paste("Best Model: ", paste(fit.min@formula)[2]))
    cat("\n------------------------------------------------------------------------------------------\n")
    return(fit.min)
  }



################################################################################


# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


# Copyrights (C) 2015, Thiago do Rego Sousa <thiagoestatistico@gmail.com>
# This is a modified version of the code contained inside file
# garch-Spec.R from package fGarch, version 3010.82.

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file



################################################################################
# FUNCTION:               SIMULATION:
#  gsSim            Simulates a GARCH/APARCH process with GEV or stable
#						        conditional distribution
################################################################################


gsSim <-
  function(spec = gsSpec(), n = 100, n.start = 100)
  {
    # A function originally implemented by Diethelm Wuertz and modified
    # to be used inside package GEVStableGarch. See the latest copyright notice.

    # Description:
    #   Simulates a time series process from the GARCH family

    # Arguments:
    #   model - a specification object of class 'GEVSTABLEGARCHSPEC' as
    #     returned by the function \code{gsSpec}:
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the intercept for ARMA specification (mean=mu/(1-sum(ar))),
    #     delta - the exponent value used in the variance
    #       equation.
    #     skew - a numeric value for the skew parameter.
    #     shape - a numeric value for the shape parameter.
    #   n - an integer, the length of the series
    #   n.start - the length of the warm-up sequence to reduce the
    #     effect of initial conditions.

    # Return:
    #   ans - An object returned by function timeDate with the simulated sample
    #     path of the specified ARMA-GARCH/APARCH model.

    # FUNCTION:

    # Error treatment of input parameters
    if(n < 2)
      stop("The parameter 'n' must be > 2")

    # Specification:
    stopifnot(class(spec) == "GEVSTABLEGARCHSPEC")
    model = spec@model

    # Random Seed:
    if (spec@rseed != 0) set.seed(spec@rseed)

    # Enlarge Series:
    n = n + n.start

    # Create Innovations:
    if (spec@distribution == "stableS0")
      z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 0)

    if (spec@distribution == "stableS1")
      z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 1)

    if (spec@distribution == "stableS2")
      z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 2)

    if (spec@distribution == "gev")
      z = rgev(n, xi = model$shape)

    if (spec@distribution == "gat")
      z = rgat(n, nu = model$shape[1], d = model$shape[2], xi = model$skew)

    if (spec@distribution == "norm")
      z = rnorm(n)

    if (spec@distribution == "std")
      z = rstd(n, nu = model$shape)

    if (spec@distribution == "sstd")
      z = rsstd(n, nu = model$shape, xi = model$skew)

    if (spec@distribution == "skstd")
      z = rskstd(n, nu = model$shape, xi = model$skew)

    if (spec@distribution == "ged")
      z = rged(n, nu = model$shape)


    # Expand to whole Sample:    NAO ENTENDI PORQUE USAR A FUNCAO rev()???
    delta = model$delta

    z = c(rev(spec@presample[, 1]), z)
    h = c(rev(spec@presample[, 2]), rep(NA, times = n))
    y = c(rev(spec@presample[, 3]), rep(NA, times = n))
    m = length(spec@presample[, 1])
    names(z) = names(h) = names(y) = NULL

    # Determine Coefficients:
    mu = model$mu
    ar = model$ar
    ma = model$ma
    omega = model$omega
    alpha = model$alpha
    gamma = model$gamma
    beta = model$beta
    deltainv = 1/delta

    # Determine Orders:
    order.ar = length(ar)
    order.ma = length(ma)
    order.alpha = length(alpha)
    order.beta = length(beta)

    # Iterate GARCH / APARCH Model and create Sample:
    # print(c(omega,alpha,gamma,beta,delta))
    eps = h^deltainv*z   # here the variable 'h' represents the process '(sigma_t)^delta'
    for (i in (m+1):(n+m)) {
      h[i] =  omega +
        sum(alpha*(abs(eps[i-(1:order.alpha)]) -
                     gamma*(eps[i-(1:order.alpha)]))^delta) +
        sum(beta*h[i-(1:order.beta)])



      eps[i] = h[i]^deltainv * z[i]
      y[i] =
        sum(ar*y[i-(1:order.ar)]) +
        sum(ma*eps[i-(1:order.ma)]) + eps[i]
    }
    y = y +  mu
    # Sample:
    data = cbind(
      z = z[(m+1):(n+m)],
      sigma = h[(m+1):(n+m)]^deltainv,
      y = y[(m+1):(n+m)])

    rownames(data) = as.character(1:n)
    if(n.start > 0)
      data = data[-(1:n.start),]


    # Return Values:
    from <-
      timeDate(format(Sys.time(), format = "%Y-%m-%d")) - NROW(data)*24*3600
    charvec  <- timeSequence(from = from, length.out = NROW(data))
    ans <- timeSeries(data = data[, c(3,2,1)], charvec = charvec)
    colnames(ans) <- c("Series", "Volatility", "Innovations")
    attr(ans, "control") <- list(gsSpec = spec)

    # Return Value:
    ans
  }



################################################################################


# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


# Copyrights (C) 2015, Thiago do Rego Sousa <thiagoestatistico@gmail.com>
# This is a modified version of the code contained inside file
# garch-Spec.R from package fGarch, version 3010.82.

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               SPECIFICATION:
#  gsSpec               Creates a 'GEVSTABLEGARCHSPEC' object from scratch
################################################################################


gsSpec <-
  function (model = list(), presample = NULL,
            cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
            rseed = NULL)
  {

    # A function originally implemented by Diethelm Wuertz and modified
    # to be used inside package GEVStableGarch. See the latest copyright notice.

    # Description:
    #   Creates a "gsSpec" object.

    # Arguments:
    #   model - a list with the model parameters as entries
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the mean value for ARMA specification,
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   presample - a numeric "matrix" with 3 columns and
    #       at least max(m,n,p,q) rows.
    #       The first culumn are the innovations, the second
    #       the conditional variances, and the last the time series.
    #       When presample is missing, we construct our presample matrix
    #       as [z,h,y] where z = rnorm(0,1), h = "uev" recursion
    #       initialization described in Wuertz et al. (2006) and
    #       y = mu. Note that the conditional variance column
    #       can contain only strictly positive numbers.
    #       If the model is a pure AR or MA the presample
    #       matrix will become a simple array.
    #   cond.dist - a character string naming the distribution
    #       function.
    #   rseed - optional random seed. The default seed is '0'.

    # Return: An object of class GEVSTABLEGARCHSPEC
    # Slots:
    #   call - the function call.
    #   formula - a formula object describing the model, e.g.
    #       ARMA(m,n) + GARCH(p,q). ARMA can be missing or
    #       specified as AR(m) or MA(n) in the case of pure
    #       autoregressive or moving average models. GARCH may
    #       alternatively be specified as ARCH(p) or APARCH(p,q).
    #   model - as declared in the input.

    # FUNCTION

    # Skewness Parameter Settings:
    skew = list(
      "stableS0" = 0,
      "stableS1" = 0,
      "stableS2" = 0,
      "gev" = NULL,
      "gat" = 1,
      "norm" = NULL,
      "std" = NULL,
      "sstd" = 0.9,
      "skstd" = 1,
      "ged" = NULL)

    # Shape Parameter Settings:
    shape = list(
      "stableS0" = 1.7,
      "stableS1" = 1.7,
      "stableS2" = 1.7,
      "gev" = 0.3,
      "gat" = c(3,1),
      "norm" = NULL,
      "std" = 4,
      "sstd" = 4,
      "skstd" = 4,
      "ged" = 2)

    # Conditional distribution
    cond.dist = match.arg(cond.dist)

    # Default Model (AR(1)):
    initialDelta = NULL
    if(!is.null(model$alpha) && is.null(model$delta)) # Garch model
    {
      if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
        initialDelta = 1
      else
        initialDelta = 2
    }

    control = list(
      omega = 1,
      alpha = NULL,
      gamma = NULL,
      beta = NULL,
      mu = NULL,
      ar = NULL,
      ma = NULL,
      delta = initialDelta,
      skew = skew[[cond.dist]],
      shape = shape[[cond.dist]]
    )

    # Update Control:
    control[names(model)] <- model
    model <- control

    # Model Orders:
    order.ar = length(model$ar)
    order.ma = length(model$ma)
    order.alpha = length(model$alpha)
    order.beta = length(model$beta)
    order.max = max(order.ar, order.ma, order.alpha, order.beta)


    # Error treatment of input parameters:
    if((order.alpha == 0 && order.beta != 0))
      stop("In a Garch(p,q)/Aparch(p,q) model we must have p > 0")

    if(!is.null(model$delta)){
      if(length(model$delta) > 1)
        stop("The parameter 'delta' must be a single number")
      if(!(model$delta > 0))
        stop("The parameter 'delta' must be > 0.")
    }

    if( (length(model$gamma) != 0) && (length(model$alpha) != 0)) # means aparch model
    {
      if(length(model$alpha) != length(model$gamma))
        stop("'alpha' and 'gamma' must have the same size for APARCH models")
    }

    if(!is.null(model$gamma)){
      if(sum(!(abs(model$gamma)<1)) > 0) # all gamma in (-1,1)
        stop("The parameter 'gamma' must be in the range -1 < gamma < 1")
    }

    if(sum(model$alpha < 0) > 0) # all alpha in [0,+infty)
      stop("The parameter 'alpha' must be greater than or equal 0")


    if(sum(model$beta < 0) > 0) # all beta in [0,+infty)
      stop("The parameter 'beta' must be greater than or equal 0")

    if(!(model$omega > 0)) # omega > 0
      stop("The parameter 'omega' must be > 0")


    if(is.null(model$ar) && is.null(model$ma) && is.null(model$alpha))
      stop("The model parameters were not specified correctly.")

    if(!is.null(model$delta) && is.null(model$alpha))
      stop("The parameter delta should be only specified for a GARCH or APARCH model.")

    # if we have a presample check if it has the correct range.
    if(is.matrix(presample) && !is.null(presample)){
      if( dim(presample)[2] != 3 || dim(presample)[1] < order.max)
        stop(cat("The presample object should be a matrix with three columns formated as: \n [Innovations, Conditional Variance, Time Series] with dimensions \n l x 3, where l = max(m,n,p,q). "))
      if( dim(presample)[1] != order.max )
      {
        warning(cat("The number of columns of the Presample matrix is \n bigger than l = max(m,n,p,q). The simulated series \n will use only the first 'l' columns"))
        presample = as.matrix ( presample[1:order.max,] )
        if(order.max == 1)
          presample = t(presample)
      }
    }

    # Compose Mean Formula Object:
    formula.mean = ""
    if (order.ar == 0 && order.ma == 0) {
      formula.mean = ""
    }
    else {
      formula.mean = paste ("arma(", as.character(order.ar), ", ",
                            as.character(order.ma), ")", sep = "")
    }

    # Compose Variance Formula Object:
    formula.var = ""
    if (order.alpha > 0) formula.var = "garch"

    if(!is.null(model$alpha)){ # decide if we have aparch model
      if (!is.null(model$gamma)){
        if(sum(model$gamma == 0) != length(model$gamma))
          formula.var = "aparch"
      } else {
        if (model$delta != 2 && ! any ( cond.dist == c("stableS0", "stableS1", "stableS2") ))
          formula.var = "aparch" # gamma = 0 and delta != 0 we get powergarch model
        if (model$delta != 1 && any ( cond.dist == c("stableS0", "stableS1", "stableS2") ))
          formula.var = "aparch"
      }
    }

    if (order.alpha == 0 && order.beta == 0) {
      formula.var = formula.var
    }
    if (order.alpha > 0) {
      formula.var = paste(formula.var, "(", as.character(order.alpha),
                          ", ", as.character(order.beta), ")", sep = "")
    }

    # Compose Mean-Variance Formula Object:
    if (formula.mean == "") {
      formula = as.formula(paste("~", formula.var))
    }
    if (formula.var == "") {
      formula = as.formula(paste("~", formula.mean))
    }
    if ((formula.mean != "") && (formula.var != "")){
      formula = as.formula(paste("~", formula.mean, "+", formula.var))
    }


    # Stop if the user specified a pure arma model
    if(formula.var == "")
      stop("Pure ARMA model not allowed")


    # Add NULL default entries:
    if (is.null(model$mu)) model$mu = 0
    if (is.null(model$ar)) model$ar = 0
    if (is.null(model$ma)) model$ma = 0
    if (is.null(model$gamma)) model$gamma = rep(0, times = order.alpha)

    # Seed:
    if (is.null(rseed)) {
      rseed = 0
    }
    else {
      set.seed(rseed)
    }

    # Define Missing Presample:
    persistency = 1-sum(model$alpha)-sum(model$beta)
    if(persistency*(1-persistency) < 0) # avoid to construct a presample with negative conditional variance.
      persistency = 0.1
    if(!is.null(presample) && is.matrix(presample)){
      z = presample[, 1]
      h = presample[, 2]
      y = presample[, 3]
      if(sum(!(h > 0)))
        stop("Conditional Variance column can have only strictly positive numbers")
    }else{
      z = rnorm(n = order.max)
      h = rep(model$omega*(1 + persistency*(1-persistency)), times = order.max)
      y = rep(model$mu, times = order.max)
    }
    presample = cbind(z, h, y)


    # Result:
    new("GEVSTABLEGARCHSPEC",
        call = match.call(),
        formula = formula,
        model = list(omega = model$omega, alpha = model$alpha,
                     gamma = model$gamma, beta = model$beta, mu = model$mu,
                     ar = model$ar, ma = model$ma, delta = model$delta,
                     skew = model$skew, shape = model$shape),
        presample = as.matrix(presample),
        distribution = as.character(cond.dist),
        rseed = as.numeric(rseed)
    )
  }


################################################################################


# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA



####################################################################################
#  FUNCTION:                        DESCRIPTION:
#
#  .stationarityAparch             Compute the stationarity condition for the
#                                   GARCH or APARCH model with several cond.
#                                   distributions
#
#  gsMomentAparch                   Evaluate the exact moments of the type
#                                   E(z - gm|z|)^delta for several conditional
#                                   distributions.
#
#  .normMomentAparch               Exact APARCH moments for the standard normal
#
#  .stdMomentAparch                Exact APARCH moments for the standard t-Student
#
#  .skstdMomentAparch              Exact APARCH moments for the standard skew
#                                   t-Student defined by Fernandez and Steel (1998).
#                                   Notice that this distribution is standardized
#                                   with location zero and unit scale, but it is not
#                                   standardized in the sense that it has a zero mean and
#                                   unit variance. The exact formula for this expression
#                                   would be very complicated to compute if we reescale
#                                   the original distribution to (X - mu) / sigma and
#                                   therefore, we choose to work with the raw distribution
#                                   defined in Fernandez and Steel without this reparametrization
#
#  .gedMomentAparch                Exact APARCH moments for the standard GED distribution

#  .gatMomentAparch                 Exact APARCH moments for the standard GAt distribution
#
#  .stableS1MomentAparch             Exact APARCH moments for the location zero and unit scale
#                                   in S1 parametrization (see Nolan (1999)).
#
#  .stableS1SymmetricMomentGarch    Exact GARCH moments ( E |zt| for the symmetric stable (S1)
#                                   distribution (S1 parametrization)
#
#  .stableS1SymmetricMomentAparch   Exact APARCH moments for the symmetric stable (S1) distribution
#
#  .stableS1MomentPowerGarch        Exact power-GARCH moments for the asymmetric stable (S1) distribution
#
#  .trueAparchMomentsWurtz           Function that evaluate moments with numerical integration.
#                                   This function was used only to test the mathematical
#                                   formulas implemented here
####################################################################################



.stationarityAparch <-
  function (model = list(),
            formula,
            cond.dist = c("stableS1", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"))
  {

    # Description:
    #   A function that evaluate stationarity conditions
    #   to guide parameter estimation in GARCH/APARCH models
    #   Make sure we have a garch or aparch with p > 0
    #   The stationarity condition for the APARCH is:
    #   sum ( E[ |z|-gm*z ] ^ delta * alpha + beta ) < 1
    #   and for the GARCH model is
    #   sum ( E[ |z| ] ^ delta * alpha + beta ) < 1
    #   where delta is usually equal to 2, except for the
    #   stable model where delta = 1.
    #   The expectation  E[ |z| ] is equal to one for the following
    #   distributions: norm, std, sstd and ged.
    #   For the GAt, gev and stable conditional distribution we need to
    #   calculate this expectation, even for the GARCH model.


    # Arguments:
    #   formula - an object returned by function .getFormula
    #   parm - a list with the model parameters as entries
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gm - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   cond.dist - a character string naming the distribution
    #       function.

    # Return Values:
    # sum ( E[ |z|-gm*z ] ^ delta * alpha + beta ) - When it can be computed
    # 1e99 - In case we the expectation E[ |z|-gm*z ] ^ delta could not be evaluated.

    # FUNCTION:

    # get parameters
    alpha <- model$alpha
    beta <- model$beta
    delta <- model$delta
    gm <- model$gm
    skew <- model$skew
    shape <- model$shape

    # Conditional distribution
    cond.dist = match.arg(cond.dist)

    #print(t(model))


    if(length(alpha) == 0  || length(alpha) != length(gm) || length(delta) != 1)
      stop("Failed to verify conditions:
           if(length(alpha) == 0 || length(delta) != 1)")

    # We must have a model with a non zero garch/aparch order
    if(formula$formula.order[3] == 0)
      stop("Invalid model")

    # garch model: sum ( E[ |z| ] ^ delta * alpha + beta ) < 1
    if(formula$isAPARCH == FALSE)
    {
      if (any( c("norm", "std", "sstd", "ged") == cond.dist))
        return(sum(alpha) + sum(beta))

      if(cond.dist == "gev")
        kappa = try(.gevMomentAparch(shape = shape, delta = 2, gm = 0), silent = TRUE)

      if(cond.dist == "stableS1")
        kappa = try(.stableS1MomentPowerGarch (shape = shape, skew = skew,
                                               delta = 1), silent = TRUE)

      if(cond.dist == "gat")
        kappa = try(.gatMomentAparch(shape = shape, skew = skew,
                                     delta = 2, gm = 0), silent = TRUE)

      if(cond.dist == "skstd")
        kappa = try(.skstdMomentAparch(shape = shape, skew = skew, delta = 2, gm = 0), silent = TRUE)

      if( is.numeric(kappa))
        return(kappa*sum(alpha) + sum(beta))
      else
        return(1e99)


    } else {
      # aparch model
      kappa = rep(0,length(alpha))

      for( i in 1:length(alpha))
      {
        if(cond.dist == "stableS1")
          kappa[i] = try(.stableS1MomentAparch (shape = shape, skew = skew,
                                                delta = delta, gm = gm[i]), silent = TRUE)
        if(cond.dist == "gev")
          kappa[i] = try(.gevMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)

        if(cond.dist == "gat")
          kappa[i] = try(.gatMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)

        if(cond.dist == "norm")
          kappa[i] = try(.normMomentAparch (delta = delta, gm = gm[i]), silent = TRUE)

        if(cond.dist == "std")
          kappa[i] = try(.stdMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)

        if(cond.dist == "skstd")
          kappa[i] = try(.skstdMomentAparch(shape = shape, skew = skew, delta = delta, gm = gm[i]), silent = TRUE)

        if(cond.dist == "ged")
          kappa[i] = try(.gedMomentAparch(shape = shape, delta = delta, gm = gm[i]), silent = TRUE)
      }
      if( is.numeric(kappa) )
        return(sum(kappa*alpha) + sum(beta))
      else
        return(1e99)
    }
  }



# ------------------------------------------------------------------------------



gsMomentAparch <- function(
    cond.dist = c("stableS1", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"),
    shape = 1.5,
    skew = 0,
    delta = 1,
    gm = 0)
{
  # conditional distribution
  cond.dist = match.arg(cond.dist)


  if(cond.dist == "stableS1")
    kappa = .stableS1MomentAparch (shape = shape, skew = skew,
                                   delta = delta, gm = gm)
  if(cond.dist == "gev")
    kappa = .gevMomentAparch(shape = shape, delta = delta, gm = gm)

  if(cond.dist == "gat")
    kappa = .gatMomentAparch(shape = shape, delta = delta, skew = skew, gm = gm)

  if(cond.dist == "norm")
    kappa = .normMomentAparch (delta = delta, gm = gm)

  if(cond.dist == "std")
    kappa = .stdMomentAparch(shape = shape, delta = delta, gm = gm)

  if(cond.dist == "sstd")
    kappa = .sstdMomentAparch(shape = shape, skew = skew, delta = delta, gm = gm)

  if(cond.dist == "skstd")
    kappa = .skstdMomentAparch(shape = shape, skew = skew, delta = delta, gm = gm)

  if(cond.dist == "ged")
    kappa = .gedMomentAparch(shape = shape, delta = delta, gm = gm)

  # Return
  kappa

}



# ------------------------------------------------------------------------------



.normMomentAparch <- function(delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard normal distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: A long memory property of stock market returns and a
  #   new model. Ding, Granger and Engle (1993), Appendix B.

  # Error treatment of input parameters
  if( ( abs ( gm ) >= 1 ) || ( delta <= 0 ) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following condition cannot be true
         ( abs ( gm ) >= 1 ) || ( delta <= 0 )")

  result = 1 / sqrt ( 2 * pi ) * ( (1 + gm ) ^ delta + ( 1 - gm ) ^ delta ) *
    2 ^ ( ( delta - 1 ) / 2 ) * gamma ( ( delta + 1 ) / 2 )

  # Return
  result
}



# ------------------------------------------------------------------------------


.stdMomentAparch <- function(shape = 3, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard t-Student distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: Mittnik -  (2000) - Conditional Density and Value-at-Risk
  #   Prediction of Asian Currency Exchange Rates
  #   Note that the formula provided by Mittnik is valid for the original t-Student
  #   distribution (which is nonStandardized).
  #   Since our objective is to calculate it for the standart t-Student distribution
  #   defined in Wurtz et al. (2006) - eq (13) we need to multiply the Mittnik
  #   formula by (sqrt((shape-2)/shape))^delta.
  #   The main drawback is that the standardized t-Student distribution
  #   is valid only for shape > 2 ( in which case the variance is finite)
  # Error treatment of input parameters
  if( (abs(gm) >= 1) || (delta <= 0) || (shape <= 2) || (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           (abs(gm) >= 1) || (delta <= 0) || (shape <= 0)|| (delta >= shape)")

  mittnikFormula = shape^(delta/2)*1/2/sqrt(pi)*( (1+gm)^delta+(1-gm)^delta )*
    gamma((delta+1)/2)*(gamma(shape/2))^(-1)*gamma((shape-delta)/2)
  factorToMultiply = (sqrt((shape-2)/shape))^delta

  # Return
  mittnikFormula*factorToMultiply
}


# ------------------------------------------------------------------------------

.skstdMomentAparch <- function(shape = 3, skew = 1, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard skew t-Student distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: fGarch: Wurtz and the package Skewt uses the one defined
  #   in Fernandez, C. and Steel, M. F. J. (1998). On Bayesian modeling of fat
  #   tails and skewness, J. Am. Statist. Assoc. 93, 359-371.
  #   Another reference: http://www.timberlake.co.uk/slaurent/G@RCH/Book63.html

  # Note: This is not the APARCH moments for the 'dsstd' function from fGarch package.
  # This is the skew t-Student defined by Fernandez and Steel without the
  # reparametrization introduced by Wurtz et al. (2006).
  # but theThe point here is that the SPLUS FinMetrics
  # also uses this distribution without reparametrization.
  # Indeed, Mittnik et al. (2000) also uses his t3-distribution (GAt) without such
  # reparametrization. These distributions are in fact standardized, but the
  # point is that the assymetry parameter plays a crucial role in the calculation
  # of the APARCH moments.

  if( (abs(gm) >= 1) || (delta <= 0) || (shape <= 2) ||
      (delta >= shape) || (skew <= 0))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           (abs(gm) >= 1) || (delta <= 0) || (shape <= 2) ||
           (delta >= shape) || (skew <= 0)")
  a = skew^(-1-delta)*(1+gm)^delta + skew^(1+delta)*(1-gm)^delta
  b = gamma((delta+1)/2)*gamma((shape-delta)/2)*(shape-2)^((1+delta)/2)
  c = (skew + 1/skew)*sqrt((shape-2)*pi)*gamma(shape/2)

  # Return
  a*b/c
}




# ------------------------------------------------------------------------------




.gatMomentAparch <- function(shape = c(3,1), skew = 1, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard GAt distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: Mittnik -  (2000) - Conditional Density and Value-at-Risk
  #   Prediction of Asian Currency Exchange Rates.

  if( (abs(gm) >= 1) || (delta <= 0) || (shape[1] <= 0) || (shape[2] <= 0) ||
      (delta >= (shape[1])*(shape[2]) ) || (skew <= 0))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true
         (abs(gm) >= 1) || (delta <= 0) || (shape[1] <= 0) || (shape[2] <= 0)
         (delta >= shape[1]*shape[2]) || (skew <= 0)")

  nu = shape[1]
  d = shape[2]
  theta = skew

  a = ( 1 + gm ) ^ delta * theta ^ ( - delta - 1 ) + ( 1 - gm ) ^ delta * theta ^ ( delta + 1 )
  b = nu ^ ( delta / d ) * beta ( ( delta + 1 ) / d , nu - delta / d )
  c = ( 1 / theta + theta ) * beta( 1 / d , nu )

  # Return
  a * b / c
}




# ------------------------------------------------------------------------------


.gedMomentAparch <- function(shape = 3, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard skew t-Student distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: http://www.timberlake.co.uk/slaurent/G@RCH/Book63.html

  if( (abs(gm) >= 1) || (delta <= 0))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
             The following conditions cannot be true
             ")
  lambda = sqrt(gamma(1/shape)*2^(-2/shape)/gamma(3/shape))
  ((1+gm)^delta + (1-gm)^delta)*2^((delta-shape)/shape)*
    gamma((delta+1)/shape)*lambda^delta/gamma(1/shape)
}



# ------------------------------------------------------------------------------



.gevMomentAparch <- function(shape = 0.3, delta = 1.2, gm = 0)
{

  # FIND A BETTER WAY TO CALCULATE IT INSTEAD OF USING THE INTEGRATION FUNCTION


  # Description:
  #   Returns the following Expectation for a standard (location zero and scale 1)
  #   GEV distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: GEVStableGarch papper

  if( ( abs( gm ) >= 1 ) || ( delta <= 0 ))
    # if( ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= -0.5 ) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           ( abs( gm ) >= 1 ) || ( delta <= 0 )")

  xi = shape

  e = function(x, gm, delta) {
    (abs(x)-gm*x)^delta * dgev(x, xi = xi)
  }

  # Compute Moment by Integration
  I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                rel.tol = .Machine$double.eps^0.25,
                gm = gm, delta = delta)

  # Return
  as.numeric(I[1])
}

# ------------------------------------------------------------------------------



.sstdMomentAparch <- function(shape = 4, skew = 1, delta = 1.2, gm = 0)
{

  # FIND A BETTER WAY TO CALCULATE IT INSTEAD OF USING THE INTEGRATION FUNCTION


  # Description:
  #   Returns the following Expectation for the sstd distribution implemented
  #   inside package fGarch
  #   sstd distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: fGarch package

  if( ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= 2 ) || ( skew <= 0 ))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
             The following conditions cannot be true
             ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= 2 ) || ( skew <= 0 )")

  nu = shape
  xi = skew

  e = function(x, gm, delta) {
    (abs(x)-gm*x)^delta * dsstd(x, nu = nu, xi = xi)
  }

  # Compute Moment by Integration
  I = integrate(e, lower = -Inf, upper = Inf, subdivisions = 1000,
                rel.tol = .Machine$double.eps^0.25,
                gm = gm, delta = delta)

  # Return
  as.numeric(I[1])
}



# ------------------------------------------------------------------------------



.stableS1MomentAparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard stable distribution in S1 parametrization
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 1.
  #   Reference: The GEVStableGarch papper on JSS.

  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) ||
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) ||
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")

  # Calculate the expression
  sigma.til <- (1 + (skew*tan(shape*pi/2))^2)^(1/2/shape)
  k.shape <- shape - 2 # since skew > 1
  skew.til <- 2/pi/(shape-2)*atan(skew*tan((shape-2)*pi/2))
  g1 <- gamma((1/2 + skew.til*k.shape/2/shape)*(-delta))
  g2 <- gamma(1/2 - skew.til*k.shape/2/shape + (1/2 + skew.til*k.shape/2/shape)*(delta+1))
  g3 <- gamma((1/2 - skew.til*k.shape/2/shape)*(-delta))
  g4 <- gamma(1/2 + skew.til*k.shape/2/shape + (1/2 - skew.til*k.shape/2/shape)*(delta+1))
  kappa <- 1/shape/sigma.til*sigma.til^(delta+1)*gamma(delta+1)*gamma(-delta/shape)*
    (1/g1/g2*(1-gm)^delta + 1/g3/g4*(1+gm)^delta)

  # Return
  kappa
}









# ------------------------------------------------------------------------------





.stableS0MomentAparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for standard stable distribution in S0 parametrization
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 0.
  #   Reference: The GEVStableGarch papper on JSS.

  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) ||
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) ||
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")

  # Compute Moment by Integration
  e = function(x, gm, delta) {
    (abs(x)-gm*x)^delta * stable::dstable.quick(x = x, alpha = shape,
                                                beta = skew, param = 0)
  }

  I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                rel.tol = .Machine$double.eps^0.25,
                gm = gm, delta = delta)

  # Return
  as.numeric(I[1])
}


.stableS2MomentAparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for standard stable distribution in S2 parametrization
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 2.
  #   Reference: The GEVStableGarch papper on JSS.

  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) ||
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) ||
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")

  # Compute Moment by Integration
  e = function(x, gm, delta) {
    (abs(x)-gm*x)^delta * stable::dstable.quick(x = x, alpha = shape,
                                                beta = skew, param = 2)
  }

  I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                rel.tol = .Machine$double.eps^0.25,
                gm = gm, delta = delta)

  # Return
  as.numeric(I[1])
}






.gevMomentAparch <- function(shape = 0.3, delta = 1.2, gm = 0)
{

  # FIND A BETTER WAY TO CALCULATE IT INSTEAD OF USING THE INTEGRATION FUNCTION


  # Description:
  #   Returns the following Expectation for a standard (location zero and scale 1)
  #   GEV distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: GEVStableGarch papper

  if( ( abs( gm ) >= 1 ) || ( delta <= 0 ))
    # if( ( abs( gm ) >= 1 ) || ( delta <= 0 ) || ( shape <= -0.5 ) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true
           ( abs( gm ) >= 1 ) || ( delta <= 0 )")

  xi = shape

  e = function(x, gm, delta) {
    (abs(x)-gm*x)^delta * dgev(x, xi = xi)
  }

  # Compute Moment by Integration
  I = integrate(e, lower = -Inf, upper = +Inf, subdivisions = 1000,
                rel.tol = .Machine$double.eps^0.25,
                gm = gm, delta = delta)

  # Return
  as.numeric(I[1])
}

# ------------------------------------------------------------------------------










# ------------------------------------------------------------------------------
.stableS1SymmetricMomentGarch <- function(shape = 1.5)
{
  # See Mittinik et al. (1995) for the definition of
  # the stable garch model with conditional symmetric stable distribution
  # This formula uses S(alpha,skew,1,0;pm) where pm = 1.
  if( (shape <= 1) || (shape > 2) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 0) || (shape > 2)")

  # Return
  gamma(1 - 1/shape)*2/pi
}



# ------------------------------------------------------------------------------


.stableS1SymmetricMomentAparch <- function(shape = 1.5, delta = 1.2, gm = 0)
{
  # See Diongue - 2008 (An investigation of the stable-Paretian Asymmetric Power GARCH model)
  # This formula uses S(alpha,0,1,0;pm) where pm = 1.
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true.
            (shape <= 1) || (shape > 2) || (delta >= shape))")


  # k(delta): At this point we have delta > 1
  if(delta == 1)
    k = pi/2
  else
    k = gamma(1-delta)*cos(pi*delta/2)
  result = k^(-1)*gamma(1-delta/shape)*1/2*((1+gm)^delta+(1-gm)^delta)

  # Return
  result
}







# ------------------------------------------------------------------------------
.stableS1MomentPowerGarch <- function(shape = 1.5, skew = 0.5, delta = 1.2)
{
  # See Mittnik et al. (2002) - Stationarity of stable power-GARCH processes
  # This formula uses S(alpha,skew,1,0;pm) where pm = 1 (tests reveald)
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) ||
      (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true.
            (shape <= 1) || (shape > 2) || ( abs(skew) > 1) ||
            (delta >= shape)")


  # k(delta): At this point we have delta > 1
  if(delta == 1)
    k = pi/2
  else
    k = gamma(1-delta)*cos(pi*delta/2)
  # tal(shape,skew)
  tau = skew*tan(shape*pi/2)
  # lambda(shape,skew,delta): eq (10) Mittnik et al. (2002)
  lambda = k^(-1)*gamma(1-delta/shape)*(1 + tau^2)^(delta/2/shape)*cos(delta/shape*atan(tau))
  lambda
}



# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
.trueAparchMomentsWurtz <-
  function(fun = "norm", gm = 0, delta = 1, lower = -Inf, upper = Inf, ...)
  {
    # This function was originally from package fGarch. We
    # changed it to return only the aparch moment expression.
    # Description:
    #   Computes APARCH moments using Integration

    # Arguments:
    #   fun - name of density functions of APARCH innovations
    #   alpha, gm - numeric value or vector of APARCH coefficients,
    #       must be of same length
    #   beta - numeric value or vector of APARCH coefficients
    #   delta - numeric value of APARCH exponent

    # Note:
    #   fun is one of: norm, snorn, std, sstd, ged, sged, snig

    # FUNCTION:

    # Match Density Function:
    fun = match.fun(fun)

    # Persisgtence Function: E(|z|-gm z)^delta
    e = function(x, gm, delta, ...) {
      (abs(x)-gm*x)^delta * fun(x, ...)
    }

    # Compute Moment by Integration
    I = integrate(e, lower = lower, upper = upper, subdivisions = 1000,
                  rel.tol = .Machine$double.eps^0.25,
                  gm = gm, delta = delta, ...)

    # Return Value:
    I
  }

################################################################################


# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


# Copyrights (C) 2015, Thiago do Rego Sousa <thiagoestatistico@gmail.com>
# This is a modified version of the code contained inside file
# garch-Spec.R from package fGarch, version 3010.82.

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file



#########################################################################################
# FUNCTION:                   DESCRIPTION:
#  show.GEVSTABLEGARCH       S4 Show method for an object of class 'GEVSTABLEGARCH'
#  show.GEVSTABLEGARCHSPEC   S4 Show method for an object of class 'GEVSTABLEGARCHSPEC'
#########################################################################################


setMethod(f = "show", signature(object = "GEVSTABLEGARCH"), definition =
            function(object)
            {
              # A function originally implemented by Diethelm Wuertz and modified
              # to be used inside package GEVStableGarch. See the latest copyright notice.

              # Description:
              #   Print method for an object of class "GEVSTABLEGARCH"

              # Arguments:
              #   object - an object of class 'GEVSTABLEGARCH'

              # FUNCTION:

              # Title:
              cat("\nTitle:\n ")
              cat(object@title, "\n")

              # Call:
              cat("\nCall:\n ")
              cat(paste(deparse(object@call), sep = "\n", collapse = "\n"), "\n")

              # Mean and Variance Equation:
              cat("\nMean and Variance Equation:\n ")
              Name = unclass(attr(object@formula, "data"))
              Formula = object@formula
              attr(Formula, "data") <- NULL
              print(Formula)
              # cat(" [", Name, "]\n", sep = "")

              # For univariate Garch Models ...

              # Conditional Distribution:
              cat("\nConditional Distribution:\n ")
              cat(object@fit$cond.dist, "\n")

              # Coefficients:
              cat("\nCoefficient(s):\n")
              digits = max(5, getOption("digits") - 4)
              print.default(format(object@fit$par, digits = digits), print.gap = 2,
                            quote = FALSE)

              # Std. Errors:
              cat("\nStd. Errors:\n ")
              cat("based on Hessian", "\n")

              # Error Analysis:
              digits = max(4, getOption("digits") - 5)
              fit = object@fit
              signif.stars = getOption("show.signif.stars")
              cat("\nError Analysis:\n")
              printCoefmat(fit$matcoef, digits = digits, signif.stars = signif.stars)

              # Log Likelihood:
              cat("\nNegative Log Likelihood:\n ")
              LLH = object@fit$llh
              N = NROW(object@data)
              cat(LLH, "   normalized: ", LLH/N, "\n")

              # Messages:
              cat("\nMessages:\n ")
              cat(object@messages$optimization.algorithm, "\n")

              # Description:
              cat("\nDescription:\n ")
              cat(object@description, "\n")

              # Return Value:
              cat("\n")
              invisible()
            })


# ------------------------------------------------------------------------------


setMethod(f = "show", signature(object = "GEVSTABLEGARCHSPEC"), definition =
            function(object)
            {
              # A function originally implemented by Diethelm Wuertz and modified
              # to be used inside package GEVStableGarch. See the latest copyright notice.

              # Description:
              #   S4 Print Method for objects of class 'GEVSTABLEGARCHSPEC'

              # Arguments:
              #   object - Object of class 'GEVSTABLEGARCHSPEC'

              # FUNCTION:

              # Formula:
              x = object
              cat("\nFormula: \n ")
              cat(as.character(x@formula))

              # Model:
              cat("\nModel:")
              if (sum(abs(x@model$ar)) != 0)
                cat("\n ar:   ", x@model$ar)
              if (sum(abs(x@model$ma)) != 0)
                cat("\n ma:   ", x@model$ma)
              if (x@model$mu != 0)
                cat("\n mu:   ", x@model$mu)
              if (x@model$omega != 0)
                cat("\n omega:", x@model$omega)
              if (sum(abs(x@model$alpha)) != 0)
                cat("\n alpha:", x@model$alpha)
              if (sum(abs(x@model$gamma)) != 0)
                cat("\n gamma:", x@model$gamma)
              if (sum(abs(x@model$beta)) != 0)
                cat("\n beta: ", x@model$beta)
              if (x@model$delta != 2)
                cat("\n delta:", x@model$delta)

              # Distribution:
              cat("\nDistribution: \n ")
              cat(x@distribution)
              if (x@distribution != "norm") {
                if (any( c("gev","std","ged") == x@distribution )) {
                  cat("\nDistributional Parameter: \n")
                  cat(" shape =", x@model$shape)
                }
                if (any( c("stable","sstd","skstd") == x@distribution )) {
                  cat("\nDistributional Parameters: \n")
                  cat(" shape =", x@model$shape, " skew =", x@model$skew)
                }
                if (x@distribution == "gat") {
                  cat("\nDistributional Parameters: \n")
                  cat(" shape = c (", x@model$shape[1], "," , x@model$shape[2], ") skew =", x@model$skew)
                }
              }

              # Seed:
              if (x@rseed != 0) {
                cat("\nRandom Seed: \n ")
                cat(x@rseed)
              }

              # Presample:
              cat("\nPresample: \n")
              n = -(length(x@presample[, 1])-1)
              time = n:0
              print(data.frame(cbind(time, x@presample)))

              # Return Value:
              invisible()
            })


################################################################################
