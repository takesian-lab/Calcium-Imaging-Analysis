function yp = resampleSINC(y,m)
%   TIME-DOMAIN SINC INTERPOLATION (RESAMPLING)
%
%   Syntax:    
%       yp = resampleSINC(y, m)
%
%   Input: 
%         y = input signal to be resampled to higher sampling rate (y must be a row vector)
%
%         m = resampling factor (e.g., if y is 100 sps and yp will be
%                                    200 sps, then m is 200/100 = 2)
%        yp = resampled signal
%
%   Example: Input is 15 Hz sinusoidal signal sampled at 200 sps. It will
%   be resampled to 400 sps using time-domain sinc interpolation
%
%       n = 256;
%       dt = 1/200;
%       t = dt*(0:n-1);
%       T = dt*n;
%       y = sin(2*pi*15*t/T);
%       m = 2;
%       yp = resampleSINC(y,m);
%       u = linspace(0,length(y),length(y));
%       up = linspace(0,length(y),length(y)*m);
%       plot(u,y,'-ob'); hold on; plot(up,yp,'-*r');
%
%   See also resampleFDZP
%    
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
%   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE
%   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
%   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
%   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 2.0 $  $Date: 2016/09/06 13:03:00 $
u = linspace(1,length(y),length(y)*m); 
x = linspace(1,length(y),length(y));
for i=1:length(u)
    yp(i) = sum(y.*sinc(u(i) - x));           
end
end
