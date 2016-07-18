mergeInto(LibraryManager.library, 
{
  databack: function(data)
  {
    self.postMessage(Module['Pointer_stringify'](data));
  },
  stamp: function()
  {
    return Math.floor(new Date().getTime() / 1000);
  }
});
