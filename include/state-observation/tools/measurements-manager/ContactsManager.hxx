namespace stateObservation
{
namespace measurements
{
///////////////////////////////////////////////////////////////////////
/// ------------------------------Contacts-----------------------------
///////////////////////////////////////////////////////////////////////

template<typename ContactT>
template<typename OnNewContact, typename OnMaintainedContact, typename OnRemovedContact, typename OnAddedContact>
void ContactsManager<ContactT>::updateContacts(const std::set<std::string> & latestContactList,
                                               OnNewContact onNewContact,
                                               OnMaintainedContact onMaintainedContact,
                                               OnRemovedContact onRemovedContact,
                                               OnAddedContact onAddedContact)
{
  contactsDetected_ = !latestContactList.empty();
  for(auto & contactName : latestContactList)
  {
    ContactT & contact = addContactToManager(contactName, onAddedContact);
    contact.wasAlreadySet(contact.isSet());
    contact.isSet(true);

    if(currentContactsList_.find(contactName) != currentContactsList_.end())
    {
      onMaintainedContact(contact);
    }
    else
    {
      onNewContact(contact);
    }
  }
  for(auto & prevContactName : currentContactsList_)
  {
    if(latestContactList.find(prevContactName) == latestContactList.end())
    {
      ContactT & removedContact = *findContact(prevContactName);
      removedContact.isSet(false);
      onRemovedContact(removedContact);
    }
  }

  currentContactsList_ = latestContactList;
}

template<typename ContactT>
template<typename OnAddedContact>
inline ContactT & ContactsManager<ContactT>::addContactToManager(const std::string & name,
                                                                 [[maybe_unused]] OnAddedContact onAddedContact)
{
  const auto [it, inserted] = listContacts_.insert({name, ContactT(idx_, name)});

  ContactT & contact = (*it).second;
  if(!inserted)
  {
    return contact;
  }

  if constexpr(!std::is_same_v<OnAddedContact, std::nullptr_t>)
  {
    onAddedContact(contact);
  }
  idx_++;

  return contact;
}
} // namespace measurements
} // namespace stateObservation
