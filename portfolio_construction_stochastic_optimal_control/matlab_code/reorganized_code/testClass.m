classdef testClass
   properties
      num;
      day;
   end
   methods
      function obj = testClass(val, day)
         obj.num = val;
         obj.day = day;
      end
      function r = print(obj)
          r = obj.num;
      end
      
      function r = getDay(obj)
          r = obj.day;
      end
      
      function r = compStr(obj, str)
          r = strcmp(str, 'CIR');
      end
   end
end

