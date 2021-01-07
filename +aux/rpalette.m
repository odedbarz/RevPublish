function RGB = rpalette(name)
% RPALETTE -- get color from Rajiv's palette.
% Usage: RGB = rpalette(name)

if isscalar(name)
    name = sprintf('new%02d', name);
end

switch (name)
   case 'black'
      RGB = [ 0 0 0 ] / 255;
      
   case 'grey'
      RGB = [ 128 128 128 ] / 255;
      
   case 'blae'
      RGB = [ 212 221 237 ] / 255;
      
   case 'white'
      RGB = [ 255 255 255 ] / 255;
      
   case 'maroon'
      RGB = [ 102 54 51 ] / 255;
      
   case 'brown'
      RGB = [ 153 102 51 ] / 255;
      
   case 'ochre'
      RGB = [ 229 164 23 ] / 255;
      
   case 'buff'
      RGB = [ 194 188 130 ] / 255;
      
   case 'yellow'
      RGB = [ 253 248 1 ] / 255;
      
   case 'scarlet'
      RGB = [ 201 45 10 ] / 255;
      
   case 'orange'
      RGB = [ 255 133 26 ] / 255;
      
   case 'peach'
      RGB = [ 255 190 140 ] / 255;
      
   case 'indigo'
      %RGB = [ 148 103 189 ] / 255;
      RGB = [ 125 45 142 ] / 255;       % Oded
      
   case 'purple'
      RGB = [ 227 119 194 ] / 255;
      
   case 'lilac'
      RGB = [ 173 143 204 ] / 255;
      
   case 'blue'
      RGB = [ 31 119 180 ] / 255;
      
   case 'azure'
      RGB = [ 117 160 191 ] / 255;
      
   case 'turquoise'
      RGB = [ 23 190 207 ] / 255;
      
   case 'teal'
      RGB = [ 6 115 84 ] / 255;
      
   case 'olive'
      RGB = [ 175 175 73 ] / 255;
      
   case 'mint'
      RGB = [ 177 235 164 ] / 255;
      
   case 'lime'
      RGB = [ 153 255 51 ] / 255;
      
   case 'forest'
      RGB = [ 44 160 44 ] / 255;
      
      
    % *** MATLAB's old lines ***
    % see: https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.bar-properties.html
    case 'line01'     
      RGB = [ 0 0 255 ] / 255;
      
    case 'line02'
      RGB = [ 0 128 0 ] / 255;
      
    case 'line03'
      RGB = [ 255 0 0 ] / 255;
      
    case 'line04'
      RGB = [ 0 191 191 ] / 255;
      
    case 'line05'
      RGB = [ 191 0 191 ] / 255;
      
    case 'line06'
      RGB = [ 191 191 0 ] / 255;
      
    case 'line07'
      RGB = [ 64 64 64 ] / 255;
      
      
    % *** MATLAB's NEW Line ***
    % see: https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.bar-properties.html
    case 'new01'    
      RGB = [0 0.4470 0.7410];
      
    case 'new02'    
      RGB = [0.8500 0.3250 0.0980];

    case 'new03'    
      RGB = [0.9290 0.6940 0.1250];
      
    case 'new04'    
      RGB = [0.4940 0.1840 0.5560];

    case 'new05'    
      RGB = [0.4660 0.6740 0.1880];
      
    case 'new06'    
      RGB = [0.3010 0.7450 0.9330   ];
      
    case 'new07'    
      RGB = [0.6350 0.0780 0.1840];
      
    otherwise
      error('Color does not belong to rpalette.');
end
