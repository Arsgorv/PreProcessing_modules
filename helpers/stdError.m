function B=stdError(A)

%try*
if sum(sum(isnan(A)))==0
    if size(A,1)==1
        B=std(A)/sqrt(length(find(~isnan(A))));
    else
        B=std(A)/sqrt(size(A,1));
    end
    
else
    if size(A,1)==1
        B=nanstd(A)/sqrt(length(find(~isnan(A))));
    else
        B=nanstd(A)/sqrt(size(A,1));
    end
end

% catch
%     B=std(A)/sqrt(length(A));
% end

