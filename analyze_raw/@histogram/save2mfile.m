function save2mfile(Obj, S, fp)


sz = get(Obj, 'Size');

Obj = struct(Obj);
Obj.size = sz;
fprintf(fp, mfileString(Obj, 'type', S));
fprintf(fp, mfileString(Obj, 'size', S));
fprintf(fp, mfileString(Obj, 'binwidth', S));
fprintf(fp, mfileString(Obj, 'offset', S));
fprintf(fp, mfileString(Obj, 'syncchan', S));
fprintf(fp, mfileString(Obj, 'spikechan', S));
fprintf(fp, mfileString(Obj, 'period', S));
fprintf(fp, mfileString(Obj, 'order', S));
fprintf(fp, mfileString(Obj, 'stimfile', S));

S = [S, '.gate'];
G = Obj.gate;
fprintf(fp, mfileString(G, 'type', S));
fprintf(fp, mfileString(G, 'delay', S, 'ms'));
fprintf(fp, mfileString(G, 'width', S, 'ms'));




