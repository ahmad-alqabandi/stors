% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proposal_info.R
\name{save_proposal}
\alias{save_proposal}
\title{Save User Proposal}
\usage{
save_proposal(proposal, proposal_name)
}
\arguments{
\item{proposal}{list representing an optimized proposal generated using the \code{build_proposal()} function.}

\item{proposal_name}{string specifying the name under which the proposal will be saved.}
}
\value{
This function will produce an error if the proposal is not generated by the \code{build_proposal()} function. Otherwise, it successfully saves the proposal without returning any value upon completion.
}
\description{
This function stores proposals generated by the \code{build_proposal()} function in R's internal data directory. It is useful when users want to reuse a proposal across multiple R sessions.
}
\examples{
# First, let's create a proposal to sample from a standard normal distribution
f_normal <- function(x) { 0.3989423 * exp( -0.5 * x^2) }
normal_proposal = build_proposal(f = f_normal, modes = 0, lower = -Inf, upper = Inf, steps = 1000)
print(normal_proposal)

# Then, we can save this proposal in R's internal data directory using `save_proposal()`
# with the name "normal"
save_proposal(normal_proposal, "normal")

# To make sure the `normal_proposal` has been stored in R's internal data directory,
# we can print all saved proposals using `print_proposals()`
print_proposals()


}
