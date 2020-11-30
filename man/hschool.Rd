\name{hschool}
\alias{hschool}
\docType{data}
\title{
Attendance Behavior of High School Juniors
}
\description{
Recorded are the number of days of absence, gender, and two test scores of 316 high school juniors from two urban high schools. 
}
\usage{data("hschool")}
\format{
  A data frame with 316 observations on the following 5 variables.
  \describe{
    \item{\code{school}}{ school of the two, coded 1 or 2;}
    \item{\code{male}}{ whether the student is male, coded 1: male; 0: female;}
    \item{\code{math}}{ the standardized test score for math;}
    \item{\code{langarts}}{ the standardized test score for language arts;}
    \item{\code{daysabs}}{ the number of days of absence.}
  }
}
\details{
Some school administrators studied the attendance behavior of high school juniors at two schools. Predictors of the number of days of absence include gender of the student and standardized test scores in math and language arts. The original source of this data set is unknown.
}
\source{
\url{https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/}
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}, \link{simglm}}.
}
\examples{
data(hschool)
summary(hschool)
head(hschool)

require(MASS)
absfit <- glm.nb(daysabs~school+male+math+langarts,data=hschool)
summary(absfit)
rsq(absfit)
rsq(absfit,adj=TRUE)

rsq.partial(absfit)
}
\keyword{datasets}
