function [feature,vector]=featureExtraction(block)

QZ=4;
FDCT=dct2(block);
m_quantized=round(FDCT./QZ);
vector=zigzag8(m_quantized);
%% reduce features
feature=vector(:,1:9);

end