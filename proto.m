obj = lfp;

ii = 1;
tcs.tt = obj.twop.ts(ii:length(obj.twop.planes.planes):end)';
dec = obj.twop.deconv(ii:length(obj.twop.planes.planes):end, obj.twop.planes.plane_members == obj.twop.planes.planes(ii));
[beh, dec] = convert_behavior(obj.behavior, tcs, dec);
analysis1 = pc_batch_analysis(beh, dec);

ii = 2;
tcs.tt = obj.twop.ts(ii:length(obj.twop.planes.planes):end)';
dec = obj.twop.deconv(ii:length(obj.twop.planes.planes):end, obj.twop.planes.plane_members == obj.twop.planes.planes(ii));
[beh, dec] = convert_behavior(obj.behavior, tcs, dec);
analysis2 = pc_batch_analysis(beh, dec);
