% axis_almost_tight()

function axis_almost_tight_positive()

axis tight
a = axis;

if sum(a(3)<0),
  a(3) = 0.001 * a(4);
end

axis([a(1)-0.05*(a(2)-a(1)), a(2)+0.05*(a(2)-a(1)),0, exp(log(a(4))+0.05*(log(a(4))-log(a(3))))])
