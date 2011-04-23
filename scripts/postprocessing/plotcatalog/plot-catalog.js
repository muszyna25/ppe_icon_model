function autoSubmit(){
  document.sidemenu.submit();
}
function findSelectedItem(selectTag)
{
for (var i=0; i<selectTag.options.length; i++)
  {
  if (selectTag.options[i].selected==true)
   {
   break
   }
  }
return i
}

function onNext(par_name)
{
var par=document.getElementById(par_name)
  i = findSelectedItem(par)
  i++
  if (i >= par.options.length)
  {
    i = 0
  }
  par.options[i].selected = true
  autoSubmit()
}

function onPrevious(par_name)
{
var par=document.getElementById(par_name)
  i = findSelectedItem(par)
  i--
  if (i < 0)
  {
    i = par.options.length - 1
  }
  par.options[i].selected = true
  autoSubmit()
}

function menuHide()
{
  document.getElementById('td_sidemenu').style.display = 'none';
  document.getElementById('menushow').style.display = 'inline';
}

function menuShow()
{
  document.getElementById('td_sidemenu').style.display = 'table-cell';
  document.getElementById('menushow').style.display = 'none';
}
