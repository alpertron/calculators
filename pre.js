var Module =
{
  'preRun': function()
  {
    self.onmessage = function(e)
	  {
  	  _doWork(Module['allocate'](Module['intArrayFromString'](e.data), 'i8', Module['ALLOC_STACK']));
	  }
  },
  'noInitialRun': true,
};

