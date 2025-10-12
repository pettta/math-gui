// boost\math\distributions\beta_binomial.hpp


// https://en.wikipedia.org/wiki/Beta-binomial_distribution
// https://www.statisticshowto.com/beta-binomial-distribution/

// The Beta-Binomial Distribution is a discrete probability distribution.
// The beta-binomial distribution is a binomial distribution in which the probability 
// of success at each trial is not fixed but randomly drawn from a beta distribution.
// It is frequently used in Bayesian statistics, empirical Bayes methods and classical 
// statistics to capture overdispersion in binomial type distributed data.
// The beta-binomial can be motivated via an urn model (Pólya urn) or as a 
// compound distribution where the success probability p ~ Beta(alpha, beta) 
// and X|p ~ Binomial(n, p).

#ifndef BOOST_MATH_DIST_BETA_BINOMIAL_HPP
#define BOOST_MATH_DIST_BETA_BINOMIAL_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/beta.hpp> // for beta function
#include <boost/math/special_functions/gamma.hpp> // for gamma function
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <boost/math/tools/roots.hpp> // for root finding

#if defined (BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4702) // unreachable code
#endif

#include <utility>

namespace boost
{
  namespace math
  {
    namespace beta_binomial_detail
    {
      // Common error checking routines for beta-binomial distribution functions:
      
      template <class RealType, class Policy>
      inline bool check_n(const char* function, const RealType& n, RealType* result, const Policy& pol)
      {
        if((n < 0) || !(boost::math::isfinite)(n))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Number of trials argument is %1%, but must be >= 0 !", n, pol);
          return false;
        }
        return true;
      } // bool check_n

      template <class RealType, class Policy>
      inline bool check_alpha(const char* function, const RealType& alpha, RealType* result, const Policy& pol)
      {
        if(!(boost::math::isfinite)(alpha) || (alpha <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Alpha argument is %1%, but must be > 0 !", alpha, pol);
          return false;
        }
        return true;
      } // bool check_alpha

      template <class RealType, class Policy>
      inline bool check_beta(const char* function, const RealType& beta, RealType* result, const Policy& pol)
      {
        if(!(boost::math::isfinite)(beta) || (beta <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Beta argument is %1%, but must be > 0 !", beta, pol);
          return false;
        }
        return true;
      } // bool check_beta

      template <class RealType, class Policy>
      inline bool check_k(const char* function, const RealType& n, const RealType& k, RealType* result, const Policy& pol)
      {
        if((k < 0) || !(boost::math::isfinite)(k))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Number of successes argument is %1%, but must be >= 0 !", k, pol);
          return false;
        }
        if(k > n)
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Number of successes argument is %1%, but must be <= Number of trials !", k, pol);
          return false;
        }
        return true;
      } // bool check_k

      template <class RealType, class Policy>
      inline bool check_prob(const char* function, const RealType& p, RealType* result, const Policy& pol)
      {
        if((p < 0) || (p > 1) || !(boost::math::isfinite)(p))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Probability argument is %1%, but must be >= 0 and <= 1 !", p, pol);
          return false;
        }
        return true;
      } // bool check_prob

      template <class RealType, class Policy>
      inline bool check_dist(const char* function, const RealType& n, const RealType& alpha, const RealType& beta, RealType* result, const Policy& pol)
      { // Check n, alpha, and beta.
        return check_n(function, n, result, pol)
          && check_alpha(function, alpha, result, pol)
          && check_beta(function, beta, result, pol);
      } // bool check_dist

      template <class RealType, class Policy>
      inline bool check_dist_and_k(const char* function, const RealType& n, const RealType& alpha, const RealType& beta, RealType k, RealType* result, const Policy& pol)
      {
        return check_dist(function, n, alpha, beta, result, pol)
          && beta_binomial_detail::check_k(function, n, k, result, pol);
      } // bool check_dist_and_k

      template <class RealType, class Policy>
      inline bool check_dist_and_prob(const char* function, const RealType& n, const RealType& alpha, const RealType& beta, RealType p, RealType* result, const Policy& pol)
      {
        return check_dist(function, n, alpha, beta, result, pol)
          && check_prob(function, p, result, pol);
      } // bool check_dist_and_prob

      template <class RealType, class Policy>
      inline bool check_mean(const char* function, const RealType& mean, RealType* result, const Policy& pol)
      {
        if(!(boost::math::isfinite)(mean) || (mean < 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "mean argument is %1%, but must be >= 0 !", mean, pol);
          return false;
        }
        return true;
      } // bool check_mean
      
      template <class RealType, class Policy>
      inline bool check_variance(const char* function, const RealType& variance, RealType* result, const Policy& pol)
      {
        if(!(boost::math::isfinite)(variance) || (variance <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "variance argument is %1%, but must be > 0 !", variance, pol);
          return false;
        }
        return true;
      } // bool check_variance
    } // namespace beta_binomial_detail

    template <class RealType = double, class Policy = policies::policy<> >
    class beta_binomial_distribution
    {
    public:
      typedef RealType value_type;
      typedef Policy policy_type;

      beta_binomial_distribution(RealType n = 1, RealType l_alpha = 1, RealType l_beta = 1) 
        : m_n(n), m_alpha(l_alpha), m_beta(l_beta)
      {
        RealType result;
        beta_binomial_detail::check_dist(
           "boost::math::beta_binomial_distribution<%1%>::beta_binomial_distribution",
          m_n,
          m_alpha,
          m_beta,
          &result, Policy());
      } // beta_binomial_distribution constructor.
      
      // Accessor functions:
      RealType trials() const
      {
        return m_n;
      }
      
      RealType alpha() const
      {
        return m_alpha;
      }
      
      RealType beta() const
      {
        return m_beta;
      }

      // Parameter estimation functions from mean and variance:
      // Based on the relationships:
      // mean = n*alpha/(alpha+beta)
      // variance = n*alpha*beta*(alpha+beta+n)/[(alpha+beta)^2*(alpha+beta+1)]

      static RealType find_alpha(
        RealType n,         // Number of trials
        RealType mean,      // Expected value of mean
        RealType variance)  // Expected value of variance
      {
        static const char* function = "boost::math::beta_binomial_distribution<%1%>::find_alpha";
        RealType result = 0;
        if(false ==
            (
              beta_binomial_detail::check_n(function, n, &result, Policy())
              && beta_binomial_detail::check_mean(function, mean, &result, Policy())
              && beta_binomial_detail::check_variance(function, variance, &result, Policy())
            )
          )
        {
          return result;
        }
        
        // From formulas: mean = n*alpha/(alpha+beta) and 
        // variance = n*p*(1-p)*[(alpha+beta+n)/(alpha+beta+1)]
        // where p = alpha/(alpha+beta)
        RealType p = mean / n;
        if(p <= 0 || p >= 1)
        {
          return policies::raise_domain_error<RealType>(
            function,
            "mean/n must be in (0,1), but is %1%", p, Policy());
        }
        
        RealType npq = n * p * (1 - p);
        RealType rho = (variance - npq) / (npq * (n - 1));
        
        if(rho <= 0)
        {
          return policies::raise_domain_error<RealType>(
            function,
            "Variance suggests underdispersion; beta-binomial not appropriate. Correlation %1%", rho, Policy());
        }
        
        RealType sum_params = (1 - rho) / rho;
        return p * sum_params;
      } // RealType find_alpha

      static RealType find_beta(
        RealType n,         // Number of trials
        RealType mean,      // Expected value of mean
        RealType variance)  // Expected value of variance
      {
        static const char* function = "boost::math::beta_binomial_distribution<%1%>::find_beta";
        RealType result = 0;
        if(false ==
            (
              beta_binomial_detail::check_n(function, n, &result, Policy())
              && beta_binomial_detail::check_mean(function, mean, &result, Policy())
              && beta_binomial_detail::check_variance(function, variance, &result, Policy())
            )
          )
        {
          return result;
        }
        
        RealType p = mean / n;
        if(p <= 0 || p >= 1)
        {
          return policies::raise_domain_error<RealType>(
            function,
            "mean/n must be in (0,1), but is %1%", p, Policy());
        }
        
        RealType npq = n * p * (1 - p);
        RealType rho = (variance - npq) / (npq * (n - 1));
        
        if(rho <= 0)
        {
          return policies::raise_domain_error<RealType>(
            function,
            "Variance suggests underdispersion; beta-binomial not appropriate. Correlation %1%", rho, Policy());
        }
        
        RealType sum_params = (1 - rho) / rho;
        return (1 - p) * sum_params;
      } // RealType find_beta

    private:
      RealType m_n;      // Number of trials
      RealType m_alpha;  // First shape parameter of beta distribution
      RealType m_beta;   // Second shape parameter of beta distribution
    }; // template <class RealType, class Policy> class beta_binomial_distribution

    #ifdef __cpp_deduction_guides
    template <class RealType>
    beta_binomial_distribution(RealType)->beta_binomial_distribution<typename boost::math::tools::promote_args<RealType>::type>;
    template <class RealType>
    beta_binomial_distribution(RealType, RealType, RealType)->beta_binomial_distribution<typename boost::math::tools::promote_args<RealType>::type>;
    #endif

    typedef beta_binomial_distribution<double> beta_binomial;

    template <class RealType, class Policy>
    inline const std::pair<RealType, RealType> range(const beta_binomial_distribution<RealType, Policy>& dist)
    { // Range of permissible values for random variable k.
      return std::pair<RealType, RealType>(static_cast<RealType>(0), dist.trials());
    }

    template <class RealType, class Policy>
    inline const std::pair<RealType, RealType> support(const beta_binomial_distribution<RealType, Policy>& dist)
    { // Range of supported values for random variable k.
      // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
      return std::pair<RealType, RealType>(static_cast<RealType>(0), dist.trials());
    }

    template <class RealType, class Policy>
    inline RealType mean(const beta_binomial_distribution<RealType, Policy>& dist)
    { // Mean of beta-binomial distribution = n*alpha/(alpha+beta).
      return dist.trials() * dist.alpha() / (dist.alpha() + dist.beta());
    } // mean

    template <class RealType, class Policy>
    inline RealType variance(const beta_binomial_distribution<RealType, Policy>& dist)
    { // Variance of beta-binomial distribution.
      // variance = n*alpha*beta*(alpha+beta+n)/[(alpha+beta)^2*(alpha+beta+1)]
      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();
      RealType sum = a + b;
      return (n * a * b * (sum + n)) / (sum * sum * (sum + 1));
    } // variance

    template <class RealType, class Policy>
    inline RealType mode(const beta_binomial_distribution<RealType, Policy>& dist)
    {
      BOOST_MATH_STD_USING // for ADL of std functions
      
      // Mode is floor((n+1)*alpha/(alpha+beta+1)) when alpha > 1
      RealType a = dist.alpha();
      RealType b = dist.beta();
      RealType n = dist.trials();
      
      if(a <= 1)
      {
        static const char* function = "boost::math::mode(beta_binomial_distribution<%1%> const&)";
        RealType result;
        result = policies::raise_domain_error<RealType>(
          function,
          "mode undefined for alpha = %1%, must be > 1!", a, Policy());
        return result;
      }
      
      // For beta-binomial: mode ≈ floor((n+1)*alpha/(alpha+beta+1))
      return floor((n + 1) * a / (a + b + 1));
    } // mode

    template <class RealType, class Policy>
    inline RealType skewness(const beta_binomial_distribution<RealType, Policy>& dist)
    {
      BOOST_MATH_STD_USING // ADL of std functions
      
      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();
      RealType sum = a + b;
      
      return ((sum + 2 * n) * (b - a) * sqrt(sum + 1)) / 
             ((sum + 2) * sqrt(n * a * b * (sum + n)));
    } // skewness

    template <class RealType, class Policy>
    inline RealType kurtosis_excess(const beta_binomial_distribution<RealType, Policy>& dist)
    {
      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();
      RealType sum = a + b;
      RealType sum2 = sum * sum;
      
      RealType numerator = (sum + 1) * (sum2 + sum * (1 - 6 * n) + 
                           n * n * (3 + 6 * n) - 
                           3 * a * b * (6 - n) / sum - 
                           18 * a * b * n * n / sum2);
      
      RealType denominator = n * a * b * (sum + 2) * (sum + 3) * (sum + n);
      
      return numerator / denominator;
    } // kurtosis_excess

    template <class RealType, class Policy>
    inline RealType kurtosis(const beta_binomial_distribution<RealType, Policy>& dist)
    {
      return 3 + kurtosis_excess(dist);
    } // kurtosis

    template <class RealType, class Policy>
    inline RealType pdf(const beta_binomial_distribution<RealType, Policy>& dist, const RealType& k)
    { // Probability Mass Function.
      BOOST_FPU_EXCEPTION_GUARD
      BOOST_MATH_STD_USING // for ADL of std functions

      static const char* function = "boost::math::pdf(beta_binomial_distribution<%1%> const&, %1%)";

      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result = 0;
      if(false == beta_binomial_detail::check_dist_and_k(
        function,
        n, a, b, k,
        &result, Policy()))
      {
        return result;
      }

      // PMF: P(X=k) = C(n,k) * B(k+alpha, n-k+beta) / B(alpha, beta)
      // Using logarithms to avoid overflow:
      // log PMF = log C(n,k) + log B(k+alpha, n-k+beta) - log B(alpha, beta)
      //         = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
      //           + lgamma(k+alpha) + lgamma(n-k+beta) - lgamma(n+alpha+beta)
      //           + lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta)
      
      using boost::math::lgamma;
      
      RealType log_result = lgamma(n + 1, Policy()) 
                          - lgamma(k + 1, Policy()) 
                          - lgamma(n - k + 1, Policy())
                          + lgamma(k + a, Policy()) 
                          + lgamma(n - k + b, Policy())
                          - lgamma(n + a + b, Policy())
                          + lgamma(a + b, Policy())
                          - lgamma(a, Policy())
                          - lgamma(b, Policy());
      
      return exp(log_result);
    } // pdf

    template <class RealType, class Policy>
    inline RealType cdf(const beta_binomial_distribution<RealType, Policy>& dist, const RealType& k)
    { // Cumulative Distribution Function beta-binomial.
      BOOST_MATH_STD_USING // for ADL of std functions

      static const char* function = "boost::math::cdf(beta_binomial_distribution<%1%> const&, %1%)";

      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result = 0;
      if(false == beta_binomial_detail::check_dist_and_k(
        function,
        n, a, b, k,
        &result, Policy()))
      {
        return result;
      }
      
      // Special cases:
      if(k == 0)
      {
        return pdf(dist, RealType(0));
      }
      if(k >= n)
      {
        return 1;
      }

      // CDF: F(k) = sum_{i=0}^{k} P(X=i)
      // We compute this as a direct sum since there's no simple closed form
      result = 0;
      for(RealType i = 0; i <= k; ++i)
      {
        result += pdf(dist, i);
      }
      
      return result;
    } // beta_binomial cdf

    template <class RealType, class Policy>
    inline RealType cdf(const complemented2_type<beta_binomial_distribution<RealType, Policy>, RealType>& c)
    { // Complemented Cumulative Distribution Function beta-binomial.
      BOOST_MATH_STD_USING // for ADL of std functions

      static const char* function = "boost::math::cdf(complement(beta_binomial_distribution<%1%> const&, %1%))";

      RealType const& k = c.param;
      beta_binomial_distribution<RealType, Policy> const& dist = c.dist;
      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result = 0;
      if(false == beta_binomial_detail::check_dist_and_k(
        function,
        n, a, b, k,
        &result, Policy()))
      {
        return result;
      }
      
      // Special cases:
      if(k >= n)
      {
        return 0;
      }
      if(k < 0)
      {
        return 1;
      }

      // Complement CDF: 1 - F(k) = sum_{i=k+1}^{n} P(X=i)
      result = 0;
      for(RealType i = k + 1; i <= n; ++i)
      {
        result += pdf(dist, i);
      }
      
      return result;
    } // beta_binomial cdf complement

    template <class RealType, class Policy>
    inline RealType quantile(const beta_binomial_distribution<RealType, Policy>& dist, const RealType& p)
    { // Quantile or Percent Point beta-binomial function.
      static const char* function = "boost::math::quantile(beta_binomial_distribution<%1%> const&, %1%)";

      RealType result = 0;
      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();
      
      if(false == beta_binomial_detail::check_dist_and_prob(
        function,
        n, a, b, p,
        &result, Policy()))
      {
        return result;
      }
      
      // Special cases:
      if(p == 0)
      {
        return 0;
      }
      if(p == 1)
      {
        return n;
      }

      // Use inverse CDF via search
      RealType cumulative = 0;
      for(RealType k = 0; k <= n; ++k)
      {
        cumulative += pdf(dist, k);
        if(cumulative >= p)
        {
          return k;
        }
      }
      
      return n;
    } // quantile

    template <class RealType, class Policy>
    inline RealType quantile(const complemented2_type<beta_binomial_distribution<RealType, Policy>, RealType>& c)
    { // Complement Quantile or Percent Point beta-binomial function.
      static const char* function = "boost::math::quantile(complement(beta_binomial_distribution<%1%> const&, %1%))";

      RealType q = c.param;
      const beta_binomial_distribution<RealType, Policy>& dist = c.dist;
      RealType result = 0;
      RealType n = dist.trials();
      RealType a = dist.alpha();
      RealType b = dist.beta();
      
      if(false == beta_binomial_detail::check_dist_and_prob(
        function,
        n, a, b, q,
        &result, Policy()))
      {
        return result;
      }
      
      // Special cases:
      if(q == 1)
      {
        return 0;
      }
      if(q == 0)
      {
        return n;
      }

      // Use inverse of complement CDF
      return quantile(dist, 1 - q);
    } // Quantile Complement

  } // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#if defined (BOOST_MSVC)
# pragma warning(pop)
#endif

#endif // BOOST_MATH_DIST_BETA_BINOMIAL_HPP
