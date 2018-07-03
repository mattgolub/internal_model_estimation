This is the codepack to accompany "Internal models for interpreting neural
population activity during sensorimotor control," by Matthew D. Golub, 
Byron M. Yu, and Steven M. Chase, eLife 2015. The original article can be 
found at: https://elifesciences.org/articles/10015/

Codepack version: 2.0
Written by: Matt Golub
Testing by: TBD

Feedback, comments, suggestions, bug reports, etc, are welcomed and 
encouraged. Please direct correspondence to mgolub@stanford.edu.

MATLAB VERSIONS:

This codepack was developed and tested using Matlab R2015a. We have also
had success with Matlab R2016a and R2017b. This codepack may not be 
compatible with earlier Matlab versions.

DESCRIPTION:

This codepack includes:
1) Three example scripts demonstrating various usages of the internal model 
estimation (IME) framework.

    example1.m:  A minimal working example of IME fitting, prediction,
                 and evaluation. This example does not use proper cross
                 validation procedures, which are critical when drawing 
                 scientific and engineering interpretations from results.
                 
                 Loads some example BMI data, fits an IME model, extracts
                 the subject's internal estimates of cursor position and 
                 velocity, and evaluate errors in the actual cursor 
                 trajectory and in the subject's internal estimates. This
                 script will also generate two figures: one quantifying 
                 those errors; and a second showing a cursor trajectory and
                 the subject's internal state estimates from an example 
                 trial. Each of these procedures can be ported for usage on 
                 new datasets as desired by leveraging 2-6 below.

    example2.m: Same as example1, but using proper cross validation 
                techniques.

    example3.m: Demonstrates model selection techniques for choosing the 
                value of TAU, the visual feedback delay used in IME models.

2) The EM algorithm for fitting an IME model.
    velime_fit.m

3) Code for extracting the subject's internal estimates of cursor
position and velocity according to the internal model.
    velime_predict.m

4) Code for evaluating and plotting the angular errors in the actual cursor
trajectories and according to the subject's internal model. Errors are 
presented in the style of Figure 3C of the paper.
    velime_evaluate.m
    plot_angular_aiming_errors.m

5) Code for visualizing single trial cursor trajectories, overlaid at
each timestep with the subject's internally predicted evolution of the 
cursor trajectory given the most recently available visual feedback of
cursor position and the subsequently issued neural commands. These
visualizations are in the style of those in Figure 4A,B of the paper.
    plot_trials_with_whiskers.m

6) Data from a representative closed-loop BMI experiment from the paper
(dataset A010509), upon which 2-5 are applied. When applying IME to new
datasets, data should be formatted according the conventions of these
example data (field names, temporal relationships, and trial structure).
    example_data.mat

AKNOWLEDGEMENTS:

This codepack leverages code written by Paolo de Leva for fast array 
expansion, multidimensional array multiplication and multidimensional array
transposing (baxfun, multiprod, and multitransp, respectively; Copyright
(c) 2009, Paolo de Leva, all rights reserved.). The license for this code 
allows for redistribution provided that the following list of conditions 
are conveyed:

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.